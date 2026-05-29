"""EnvironmentManager.check_environments() caching (audit 2026-05-28,
F-ENG-PERF). Enumerating conda envs spawns a subprocess; run_in_stack used to
do it on every dispatch. The result must be cached and reused while the
env-aliases file is unchanged."""

import json


def test_check_environments_is_cached(tmp_path, monkeypatch):
    from aria.utils import environment_manager as em

    mgr = em.EnvironmentManager.__new__(em.EnvironmentManager)
    mgr._conda_ok = True
    mgr._env_cache = None
    mgr._env_cache_key = None

    calls = {"n": 0}

    class _Result:
        stdout = json.dumps({"envs": ["/x/envs/aria-rna-env"]})

    def fake_run(cmd, *a, **k):
        if cmd[:3] == ["conda", "env", "list"]:
            calls["n"] += 1
        return _Result()

    monkeypatch.setattr(em.subprocess, "run", fake_run)
    # Force a stable cache key (no aliases file) by pointing HOME at an empty dir.
    monkeypatch.setenv("HOME", str(tmp_path))

    first = mgr.check_environments()
    second = mgr.check_environments()

    assert first == second
    assert first["rna"] is True
    assert calls["n"] == 1, "second call should hit the cache, not subprocess"
