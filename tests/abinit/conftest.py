import logging
from pathlib import Path
from typing import Literal, Sequence, Union

import pytest

logger = logging.getLogger("atomate2")

_REF_PATHS = {}
_ABINIT_FILES = (
    "run.abi",
    "abinit_input.json",
)
_FAKE_RUN_ABINIT_KWARGS = {}


@pytest.fixture(scope="session")
def abinit_test_dir(test_dir):
    return test_dir / "abinit"


@pytest.fixture(scope="session")
def abinit_integration_tests(pytestconfig):
    return pytestconfig.getoption("abinit_integration")


@pytest.fixture(scope="function")
def mock_abinit(mocker, abinit_test_dir, abinit_integration_tests):
    """
    This fixture allows one to mock running ABINIT.

    It works by monkeypatching (replacing) calls to run_abinit.

    The primary idea is that instead of running ABINIT to generate the output files,
    reference files will be copied into the directory instead.
    """
    import atomate2.abinit.files
    import atomate2.abinit.jobs.base
    import atomate2.abinit.run

    # Wrap the write_abinit_input_set so that we can check inputs after calling it
    def wrapped_write_abinit_input_set(*args, **kwargs):
        from jobflow import CURRENT_JOB

        name = CURRENT_JOB.job.name
        ref_path = abinit_test_dir / _REF_PATHS[name]

        atomate2.abinit.files.write_abinit_input_set(*args, **kwargs)
        check_abinit_inputs(ref_path)

    mocker.patch.object(
        atomate2.abinit.jobs.base,
        "write_abinit_input_set",
        wrapped_write_abinit_input_set,
    )

    if not abinit_integration_tests:
        # Mock abinit run (i.e. this will copy reference files)
        def mock_run_abinit(wall_time=None, start_time=None):
            from jobflow import CURRENT_JOB

            name = CURRENT_JOB.job.name
            ref_path = abinit_test_dir / _REF_PATHS[name]
            check_abinit_inputs(ref_path)
            fake_run_abinit(ref_path)

        mocker.patch.object(atomate2.abinit.run, "run_abinit", mock_run_abinit)
        mocker.patch.object(atomate2.abinit.jobs.base, "run_abinit", mock_run_abinit)

        def _run(ref_paths, fake_run_abinit_kwargs=None):
            if fake_run_abinit_kwargs is None:
                fake_run_abinit_kwargs = {}
            _REF_PATHS.update(ref_paths)
            _FAKE_RUN_ABINIT_KWARGS.update(fake_run_abinit_kwargs)

        yield _run

    mocker.stopall()
    _REF_PATHS.clear()
    _FAKE_RUN_ABINIT_KWARGS.clear()


def fake_run_abinit(
    ref_path: Union[str, Path],
):
    """
    Emulate running ABINIT.

    Parameters
    ----------
    ref_path
        Path to reference directory with ABINIT input files in the folder named 'inputs'
        and output files in the folder named 'outputs'.
    """
    logger.info("Running fake ABINIT.")

    ref_path = Path(ref_path)

    copy_abinit_outputs(ref_path)

    # pretend to run ABINIT by copying pre-generated outputs from reference dir
    logger.info("Generated fake ABINIT outputs")


def check_abinit_inputs(
    ref_path: Union[str, Path],
    check_inputs: Sequence[Literal["run.abi"]] = _ABINIT_FILES,
):
    ref_path = Path(ref_path)

    if "run.abi" in check_inputs:
        check_run_abi(ref_path)

    if "abinit_input.json" in check_inputs:
        check_abinit_input_json(ref_path)

    logger.info("Verified inputs successfully")


def check_run_abi(ref_path: Union[str, Path]):
    from abipy.abio.abivars import AbinitInputFile

    user = AbinitInputFile.from_file("run.abi")
    assert user.ndtset == 1, f"'run.abi' has multiple datasets (ndtset={user.ndtset})."
    ref = AbinitInputFile.from_file(ref_path / "inputs" / "run.abi")
    # Ignore the pseudos as the directory depends on the pseudo root directory
    diffs = user.get_differences(ref, ignore_vars=["pseudos"])
    # TODO: should we still add some check on the pseudos here ?
    assert diffs == [], "'run.abi' is different from reference."


def check_abinit_input_json(ref_path: Union[str, Path]):
    from abipy.abio.inputs import AbinitInput
    from monty.serialization import loadfn

    user = loadfn("abinit_input.json")
    assert isinstance(user, AbinitInput)
    ref = loadfn(ref_path / "inputs" / "abinit_input.json")
    assert user.structure == ref.structure
    assert user.runlevel == ref.runlevel


def clear_abinit_files():
    for abinit_file in ("run.abo",):
        if Path(abinit_file).exists():
            Path(abinit_file).unlink()
    logger.info("Cleared abinit files.")


def copy_abinit_outputs(ref_path: Union[str, Path]):
    import shutil

    output_path = ref_path / "outputs"
    for output_file in output_path.iterdir():
        if output_file.is_file():
            shutil.copy(output_file, ".")
