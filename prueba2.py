import time, sys, uuid
sys.path.insert(0, '.')

from aria.memory.memory import ARIAMemory
from aria.llm.provider import LLMProvider
from aria.agents.bulk_rna_agent import BulkRNAAgent, _is_fastq
from aria.bus.message_bus import bus

memory = ARIAMemory(':memory:')
llm    = LLMProvider.from_config()
agent  = BulkRNAAgent(memory, llm)
exp_id = 'test_' + uuid.uuid4().hex[:6]

files = [
    '/home/medusa/Samael/H9-RNA/raw_fastq/B1_1.fq.gz',
    '/home/medusa/Samael/H9-RNA/raw_fastq/B1_2.fq.gz',
    '/home/medusa/Samael/H9-RNA/raw_fastq/WT1_1.fq.gz',
    '/home/medusa/Samael/H9-RNA/raw_fastq/WT1_2.fq.gz',
]

print('is_fastq:', _is_fastq(files))
print('Starting run...')

import threading
done = threading.Event()

def run():
    result = agent.run(exp_id, {
        'exp_context': {
            'organism':  'Homo sapiens',
            'genome':    'hg38',
            'modalities': {'bulk_RNA': files},
            'genome_config': {
                'genome_key':  'hg38',
                'fasta':       str(__import__('pathlib').Path.home() / '.aria/genomes/hg38/genome.fa.gz'),
                'gtf':         str(__import__('pathlib').Path.home() / '.aria/genomes/hg38/annotation.gtf.gz'),
                'star_index':  str(__import__('pathlib').Path.home() / '.aria/genomes/hg38/star_index'),
                'strand': 0,
            }
        },
        'biological_intent': {
            'comparison': 'BMAL1 KO vs WT',
            'summary':    'DE genes BMAL1 KO vs WT',
        }
    })
    print('run() returned:', result.get('status'))
    print('keys:', list(result.keys()))
    done.set()

t = threading.Thread(target=run, daemon=True)
t.start()

# Print progress every 3 seconds
for i in range(20):
    if done.wait(timeout=3):
        break
    print(f'  ...still running ({i*3}s elapsed)')
    findings = bus.get_findings(exp_id)
    if findings:
        print(f'  findings so far: {len(findings)}')

print('done')

