"""Read-only checks of saved strengthening analyses, without fitting any model.

The separate synthetic two-tree scaler-guard unit test is deliberately omitted;
all analyzed fitted models, their training-only scaling, selections, predictions
and independent metrics are still checked. No report or bytecode is written.
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError('Run without Python optimization; verification assertions are mandatory')
sys.dont_write_bytecode=True

import model_verify
import verify_strengthening as independent
import verify_mlp_budget_sensitivity as sensitivity
from prediction_portability_policy import validate_portable_diagnostic


def verify(portable_float32=False,artifact_manifest=None):
    if artifact_manifest is not None and not portable_float32:raise ValueError('An alternate diagnostic manifest requires explicit portable inference')
    # This wrapper calls only the no-write interfaces, never their report-writing
    # command-line entry points. Captured test output is emitted only on failure.
    captured=io.StringIO()
    try:
        with contextlib.redirect_stdout(captured),contextlib.redirect_stderr(captured):
            models=model_verify.main(write=False,include_toy_fit=False,portable_float32=portable_float32)
        assert models['passed'] and models['artifact_reload_predictions']==164
        assert models['verification_files_written']==0 and not models['toy_scaler_guard_fit_executed']
        checks={
            'models':models,
            'doublets':independent.audit_doublets(),
            'second_source':independent.audit_models(independent.DATA/'models/secondsource',portable_float32=portable_float32),
            'nested_model_comparison':independent.audit_models(independent.DATA/'models/nested_comparison',nested=True,portable_float32=portable_float32),
            'mlp_budget_sensitivity':sensitivity.verify(write=False,portable_float32=portable_float32),
        }
        assert all(checks[k]['complete'] for k in ['doublets','second_source','nested_model_comparison'])
        assert checks['mlp_budget_sensitivity']['passed']
        assert checks['mlp_budget_sensitivity']['outer_models_reloaded']==10
        if portable_float32:
            from diagnose_prediction_portability import diagnose
            with contextlib.redirect_stdout(captured),contextlib.redirect_stderr(captured):
                diagnostic=diagnose(manifest=artifact_manifest)
            checks['portable_float32_inference']=validate_portable_diagnostic(diagnostic)
        checks.update(passed=True,model_fits_started=0,files_written=0,
            scope='Saved scientific analyses and existing fitted models; no report regeneration or synthetic guard fit')
        return checks
    except Exception:
        sys.stderr.write(captured.getvalue())
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--portable-float32',action='store_true',help='Explicit cross-platform float32 MLP inference check with mandatory metric/ranking/decision guards; no archived score changes')
    parser.add_argument('--artifact-manifest',type=Path,help='Local diagnostic input inventory only; public extracted snapshots already contain their exact manifest')
    args=parser.parse_args()
    print(json.dumps(verify(portable_float32=args.portable_float32,artifact_manifest=args.artifact_manifest),indent=2))
