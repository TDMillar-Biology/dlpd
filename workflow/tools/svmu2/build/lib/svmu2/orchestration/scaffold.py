'''
Docstring for orchestration.scaffold
Accept args from cli layer to run scaffolding logic
'''
from svmu2.orchestration.parse import run_parse
from svmu2.core.scaffold import atomicity_filter

def run_scaffold(args):
    all_alns, primary = run_parse(args)
    atomicity_filter(all_alns)