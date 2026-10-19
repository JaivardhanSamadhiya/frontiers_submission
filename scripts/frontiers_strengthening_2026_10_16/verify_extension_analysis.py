"""Read-only checks of saved strengthening analyses, without fitting any model.

The separate synthetic two-tree scaler-guard unit test is deliberately omitted;
all analyzed fitted models, their training-only scaling, selections, predictions
and independent metrics are still checked. No report or bytecode is written.
"""
from __future__ import annotations

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


def verify():
    # This wrapper calls only the no-write interfaces, never their report-writing
    # command-line entry points. Captured test output is emitted only on failure.
    captured=io.StringIO()
    try:
        with contextlib.redirect_stdout(captured),contextlib.redirect_stderr(captured):
            models=model_verify.main(write=False,include_toy_fit=False)
        assert models['passed'] and models['artifact_reload_predictions']==164
        assert models['verification_files_written']==0 and not models['toy_scaler_guard_fit_executed']
        checks={
            'models':models,
            'doublets':independent.audit_doublets(),
            'second_source':independent.audit_models(independent.DATA/'models/secondsource'),
            'nested_model_comparison':independent.audit_models(independent.DATA/'models/nested_comparison',nested=True),
            'mlp_budget_sensitivity':sensitivity.verify(write=False),
        }
        assert all(checks[k]['complete'] for k in ['doublets','second_source','nested_model_comparison'])
        assert checks['mlp_budget_sensitivity']['passed']
        assert checks['mlp_budget_sensitivity']['outer_models_reloaded']==10
        checks.update(passed=True,model_fits_started=0,files_written=0,
            scope='Saved scientific analyses and existing fitted models; no report regeneration or synthetic guard fit')
        return checks
    except Exception:
        sys.stderr.write(captured.getvalue())
        raise


if __name__=='__main__':print(json.dumps(verify(),indent=2))
