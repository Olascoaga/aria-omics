"""Typed verification receipt for the SPECULATIVE causal gate (ADR-057 rail #1).

Round-3 H14 (Codex blocker 1). ADR-057 rail #1 says the agent only speculates
downstream of W-CLAIM + W-LEDGER PASSING. The previous design carried that state
as two bare booleans synthesised by ``report_builder._speculative_verification_state``
— and on ABSENCE of the verification artifacts that function returned
``(True, True)``, i.e. *absence was read as approval*. That made the wall
decorative exactly in the dangerous case (a run whose verification did not
complete would still get speculation layered on top).

``VerificationReceipt`` replaces the boolean pair with a typed, fail-closed
value. It carries not only the two pass flags but ``complete`` — did we actually
OBSERVE both verification artifacts? — plus a short human-readable evidence
string for each side. The gate opens ONLY when the receipt is complete AND both
sides passed. Absence ⇒ ``complete=False`` ⇒ gate closed, with a distinct,
auditable reason (``verification_evidence_absent``) so an honest-null is never
confused with a silent fail-open.

Pure data: asserts no biology, runs no LLM.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

# Reasons the gate stayed shut (surfaced in the agent result + manifest).
VERIFICATION_EVIDENCE_ABSENT = "verification_evidence_absent"
VERIFICATION_GATE_NOT_PASSED = "verification_gate_not_passed"


@dataclass(frozen=True)
class VerificationReceipt:
    """Fail-closed record of the run's W-CLAIM / W-LEDGER state.

    ``complete`` is the fail-closed pivot: it is True only when BOTH verification
    artifacts were actually observed. A receipt built from absence is incomplete
    and never opens the gate, regardless of the (defaulted) pass flags.
    """

    w_claim_passed: bool
    w_ledger_passed: bool
    complete: bool
    w_claim_evidence: str = ""
    w_ledger_evidence: str = ""

    @property
    def gate_open(self) -> bool:
        """True iff verification is complete AND both sides passed."""
        return bool(self.complete and self.w_claim_passed and self.w_ledger_passed)

    @property
    def blocked_reason(self) -> str:
        """Why the gate stayed shut ('' when it is open).

        Distinguishes *we never saw the verification* (evidence absent — the
        round-3 blocker) from *we saw it and it failed* (gate not passed).
        """
        if self.gate_open:
            return ""
        if not self.complete:
            return VERIFICATION_EVIDENCE_ABSENT
        return VERIFICATION_GATE_NOT_PASSED

    def to_dict(self) -> dict:
        data = asdict(self)
        data["gate_open"] = self.gate_open
        data["blocked_reason"] = self.blocked_reason
        return data

    @classmethod
    def from_explicit(
        cls, w_claim_passed: bool, w_ledger_passed: bool
    ) -> "VerificationReceipt":
        """Build a complete receipt from a caller's explicit assertion.

        An explicit ``w_claim_passed`` / ``w_ledger_passed`` IS positive evidence:
        the caller (or a test) is asserting the run's resolved state, so the
        receipt is ``complete``. This is the backward-compatible path for callers
        that already resolved verification to two booleans (ADR-057 H13). The
        fail-open hazard lived in *synthesising* a pass from absence, which this
        path never does — a caller that asserts nothing gets no receipt at all.
        """
        return cls(
            w_claim_passed=bool(w_claim_passed),
            w_ledger_passed=bool(w_ledger_passed),
            complete=True,
            w_claim_evidence="explicit caller assertion",
            w_ledger_evidence="explicit caller assertion",
        )
