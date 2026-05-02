#!/usr/bin/env bash
# ============================================================
# ARIA v4.0 — Fix: paired-end merging, contrast injection
# ============================================================
set -euo pipefail

REPO_DIR="/home/medusa/Samael/ARIA"
BACKUP_DIR="/tmp/aria_backup_$(date +%s)"
mkdir -p "$BACKUP_DIR"
mkdir -p "$BACKUP_DIR/aria/agents"

echo "==> Backing up files to $BACKUP_DIR"
for f in aria/agents/design_agent.py aria/agents/orchestrator_agent.py aria/agents/bulk_rna_agent.py; do
    cp -v "$REPO_DIR/$f" "$BACKUP_DIR/$f"
done

# ── 1. design_agent.py: _parse_samples paired-end merging ──
# Reemplaza la función existente con la nueva versión.
# La nueva función agrupa reads R1/R2.
echo "==> Patching design_agent.py (_parse_samples)"

awk '
/def _parse_samples\(self, file_paths: list\[str\]\) -> list\[dict\]:/ {
    in_func = 1
    in_indent = 0
    print
    print "        \"\"\"Parse sample names, merging paired-end read pairs (R1/R2).\"\"\""
    print "        # Agrupar por nombre base (quitando _1/_2 o _R1/_R2)"
    print "        pairs = {}"
    print "        for f in file_paths:"
    print "            stem = Path(f).stem"
    print "            while stem.endswith((\".fastq\", \".fq\", \".gz\", \".tsv\", \".csv\", \".txt\")):"
    print "                stem = Path(stem).stem"
    print "            # Normalizar nombres paired-end: quitar sufijo _1/_2 o _R1/_R2"
    print "            base = re.sub(r\"(_[12]|_R[12])$\", \"\", stem, flags=re.IGNORECASE)"
    print "            pairs.setdefault(base, []).append(f)"
    print ""
    print "        parsed = []"
    print "        for base, files in pairs.items():"
    print "            tokens = re.split(r\"[_\\-\\.\\s]+\", base.lower())"
    print "            tokens = [t for t in tokens if t]"
    print "            parsed.append({"
    print "                \"raw\": files[0],         # usar el primer archivo del par"
    print "                \"stem\": base,"
    print "                \"tokens\": tokens,"
    print "                \"paired_files\": files    # todos los archivos del par"
    print "            })"
    print "        return parsed"
    next
}
in_func && /^    def / { in_func=0 }
in_func && /^    @staticmethod/ { print; next }
in_func { next }
{ print }
' "$REPO_DIR/aria/agents/design_agent.py" > /tmp/design_agent_patched.py

mv /tmp/design_agent_patched.py "$REPO_DIR/aria/agents/design_agent.py"

# ── 2. orchestrator_agent.py: solicitar contrastes en plan, inyectarlos en diseño ──
echo "==> Patching orchestrator_agent.py (plan contrasts)"

# 2a. Modificar _design_analysis_plan para incluir campo "contrasts"
awk '
/def _design_analysis_plan\(self, experiment_id: str, exp_context: dict\) -> dict:/ {
    print
    print "        prompt = f\"\"\""
    print "Available modalities: {list(exp_context.get(\"modalities\", {}).keys())}"
    print "User question: \"{exp_context.get(\"user_question\", \"\")}\""
    print "Organism: {exp_context.get(\"organism\")}"
    print "Genome: {exp_context.get(\"genome\")}"
    print "Is multimodal: {exp_context.get(\"is_multimodal\")}"
    print ""
    print "Design the analysis pipeline. Return JSON:"
    print "{{"
    print "  \"steps\": [{{\"order\": 1, \"agent\": \"rna_agent\", \"analysis\": \"description\", \"depends_on\": [], \"can_parallel\": false}}],"
    print "  \"contrasts\": [{{\"numerator\": \"groupA\", \"denominator\": \"groupB\"}}],"
    print "  \"integration_needed\": true,"
    print "  \"integration_type\": \"WNN|MOFA|none\","
    print "  \"estimated_complexity\": \"low|medium|high\","
    print "  \"rationale\": \"brief explanation\""
    print "}}"
    print "\"\"\""
    in_func=1
    next
}
in_func && /^    def / { in_func=0 }
in_func { next }
{ print }
' "$REPO_DIR/aria/agents/orchestrator_agent.py" > /tmp/orch1.py

# 2b. Inyectar contrastes del plan en el diseño (dentro de _after_checkpoint_2)
awk '
/        self\.memory\.store_decision\(/ && /Analysis plan confirmed/ {
    print
    print ""
    print "        # v4.0: inject plan contrasts into design for BulkRNAAgent"
    print "        plan_contrasts = plan.get(\"contrasts\")"
    print "        if plan_contrasts and exp_context.get(\"design\"):"
    print "            exp_context[\"design\"][\"plan_contrasts\"] = plan_contrasts"
    print ""
    print "        threading.Thread("
    in_inject=1
    next
}
in_inject && /threading.Thread/ { in_inject=0 }
{ print }
' /tmp/orch1.py > /tmp/orch2.py

mv /tmp/orch2.py "$REPO_DIR/aria/agents/orchestrator_agent.py"

# ── 3. bulk_rna_agent.py: mejorar _apply_design y fusionar contrastes del plan ──
echo "==> Patching bulk_rna_agent.py (apply design + plan contrasts)"

# 3a. Dentro de _apply_design, antes del fallback "if not group_labels:", intentar recorte de sufijos
awk '
/        # Fallback: if no match, assign from design directly \(assumes column order\)/ {
    print
    print "        # Attempt to match by trimming potential technical replicate suffixes (_1, _2, _R1, _R2)"
    print "        if not group_labels:"
    print "            for col in ordered_cols:"
    print "                best_match = None"
    print "                for stem, grp in sample_stems.items():"
    print "                    # Trim common suffixes from col"
    print "                    col_base = re.sub(r\"[_\\-]?[12]$\", \"\", col)"
    print "                    if col_base == stem or stem.startswith(col_base) or col_base.startswith(stem):"
    print "                        if best_match is not None:"
    print "                            best_match = None"
    print "                            break"
    print "                        best_match = grp"
    print "                if best_match:"
    print "                    group_labels[col] = best_match"
    in_fallback=1
    next
}
in_fallback && /        if not group_labels or len/ { in_fallback=0 }
{ print }
' "$REPO_DIR/aria/agents/bulk_rna_agent.py" > /tmp/bulk1.py

# 3b. Después de construir contrastes (tanto design como legacy), añadir los del plan
awk '
/        self\.publish_status\(/ && /Detected groups/ {
    print
    print "        # Add plan contrasts if available"
    print "        if design and design.get(\"plan_contrasts\"):"
    print "            existing = set((c[\"numerator\"], c[\"denominator\"]) for c in contrasts)"
    print "            for pc in design[\"plan_contrasts\"]:"
    print "                num, den = pc[\"numerator\"], pc[\"denominator\"]"
    print "                if (num, den) not in existing and (den, num) not in existing:"
    print "                    contrasts.append(pc)"
    print "                    existing.add((num, den))"
    print ""
    in_contrasts=1
    next
}
in_contrasts && /^        self\.publish_status/ { in_contrasts=0 }
{ print }
' /tmp/bulk1.py > /tmp/bulk2.py

mv /tmp/bulk2.py "$REPO_DIR/aria/agents/bulk_rna_agent.py"

echo "==> Done. All patches applied."
echo "    Backup at: $BACKUP_DIR"
