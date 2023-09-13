import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import pytest
from abipy.abio.input_tags import SCF
from abipy.abio.inputs import AbinitInput
from monty.os import makedirs_p
from monty.tempfile import ScratchDir

import atomate2
from atomate2.abinit.files import load_abinit_input
from atomate2.abinit.sets.base import (
    AbinitInputGenerator,
    AbinitInputSet,
    as_pseudo_table,
)
from atomate2.abinit.utils.common import INDIR_NAME, OUTDIR_NAME, InitializationError


class TestAbinitInputSet:
    def test_init(self, abinit_test_dir):
        abinit_input = load_abinit_input(
            os.path.join(abinit_test_dir, "abinit_inputs"), fname="abinit_input_Si.json"
        )
        ais = AbinitInputSet(abinit_input=abinit_input)
        assert set(ais.inputs.keys()) == {"run.abi", "abinit_input.json"}
        assert ais.inputs["run.abi"] == abinit_input
        assert '"@class": "AbinitInput"' in ais.inputs["abinit_input.json"]
        assert ais.input_files is None
        assert ais.link_files is True
        assert ais.validate() is True
        ais = AbinitInputSet(
            abinit_input=abinit_input,
            input_files=[
                ("/some/input/file", "file"),
                ("/some/other/input/file", "other_file"),
            ],
            link_files=False,
        )
        assert ais.input_files == [
            ("/some/input/file", "file"),
            ("/some/other/input/file", "other_file"),
        ]
        assert ais.link_files is False
        assert ais.validate() is False

    def test_set_workdir(self):
        with ScratchDir(".") as tmp:
            indir, outdir, tmpdir = AbinitInputSet.set_workdir("someworkdir")
            indata = os.path.join("someworkdir", "indata")
            outdata = os.path.join("someworkdir", "outdata")
            tmpdata = os.path.join("someworkdir", "tmpdata")
            assert os.path.exists(indata)
            assert os.path.exists(outdata)
            assert os.path.exists(tmpdata)
            assert indir.path == os.path.join(tmp, indata)
            assert outdir.path == os.path.join(tmp, outdata)
            assert tmpdir.path == os.path.join(tmp, tmpdata)

    def test_write_input(self, abinit_test_dir):
        abinit_input = load_abinit_input(
            os.path.join(abinit_test_dir, "abinit_inputs"), fname="abinit_input_Si.json"
        )
        with ScratchDir("."):
            ais = AbinitInputSet(abinit_input=abinit_input)
            assert ais.validate() is True
            ais.write_input("testdir")
            assert os.path.exists("testdir")
            dirlist = os.listdir("testdir")
            assert len(dirlist) == 5
            assert "indata" in dirlist
            assert "outdata" in dirlist
            assert "tmpdata" in dirlist
            assert os.path.isdir("testdir/indata")
            assert os.path.isdir("testdir/outdata")
            assert os.path.isdir("testdir/tmpdata")
            assert "run.abi" in dirlist
            assert "abinit_input.json" in dirlist
            with open("testdir/run.abi", "r") as f:
                abistr = f.read()
                assert "ecut" in abistr
            with open("testdir/abinit_input.json", "r") as f:
                abijsonstr = f.read()
                assert "@module" in abijsonstr
            with pytest.raises(FileExistsError):
                ais.write_input("testdir", overwrite=False)

        with ScratchDir(".") as tmpdir:
            prev_output_dir = os.path.join(tmpdir, "prev_output")
            prev_outdata = os.path.join(prev_output_dir, OUTDIR_NAME)
            makedirs_p(prev_outdata)
            out_den = Path(os.path.join(prev_outdata, "out_DEN"))
            out_den.touch()
            ais = AbinitInputSet(
                abinit_input=abinit_input, input_files=[(out_den, "in_DEN")]
            )
            assert "irdden" not in ais.abinit_input
            assert ais.validate() is False
            ais.write_input("testdir")
            in_den = os.path.join("testdir", "indata", "in_DEN")
            assert os.path.islink(in_den)
            assert os.readlink(in_den) == str(out_den)

        with ScratchDir(".") as tmpdir:
            prev_output_dir = os.path.join(tmpdir, "prev_output")
            prev_outdata = os.path.join(prev_output_dir, OUTDIR_NAME)
            makedirs_p(prev_outdata)
            out_den = Path(os.path.join(prev_outdata, "out_DEN"))
            out_den.touch()
            ais = AbinitInputSet(
                abinit_input=abinit_input,
                input_files=[(out_den, "in_DEN")],
                link_files=False,
            )
            assert "irdden" not in ais.abinit_input
            assert ais.validate() is False
            ais.abinit_input["irdden"] = 1
            assert "irdden" in ais.abinit_input
            assert ais.validate() is True
            ais.write_input("testdir")
            assert os.path.exists(in_den)
            assert os.path.isfile(in_den)
            assert not os.path.islink(in_den)
            with open("testdir/run.abi", "r") as f:
                abistr = f.read()
                assert "irdden 1" in abistr
            del ais.abinit_input["irdden"]

        with ScratchDir(".") as tmpdir:
            prev_output_dir = os.path.join(tmpdir, "prev_output")
            prev_outdata = os.path.join(prev_output_dir, OUTDIR_NAME)
            makedirs_p(prev_outdata)
            Path(os.path.join(prev_outdata, "out_TIM1_DEN")).touch()
            Path(os.path.join(prev_outdata, "out_TIM2_DEN")).touch()
            out_den = Path(os.path.join(prev_outdata, "out_TIM15_DEN"))
            out_den.touch()
            abinit_input["irdden"] = 1
            ais = AbinitInputSet(
                abinit_input=abinit_input,
                input_files=[(str(out_den), "in_DEN")],
            )
            ais.write_input("testdir")
            in_den = os.path.join("testdir", INDIR_NAME, "in_DEN")
            assert os.path.islink(in_den)
            assert os.readlink(in_den) == str(out_den)

    def test_set_remove_vars(self, abinit_test_dir):
        abinit_input = load_abinit_input(
            os.path.join(abinit_test_dir, "abinit_inputs"), fname="abinit_input_Si.json"
        )
        ais = AbinitInputSet(abinit_input=abinit_input)
        assert "ntime" not in ais.abinit_input
        assert ais.abinit_input["nband"] == 5
        ais.set_vars(nband=10, ntime=5)
        assert ais.abinit_input["nband"] == 10
        assert ais.abinit_input["ntime"] == 5
        assert "occopt" in ais.abinit_input
        ais.remove_vars("occopt")
        assert "occopt" not in ais.abinit_input
        ais.remove_vars(["nband", "ntime"])
        assert "nband" not in ais.abinit_input
        assert "ntime" not in ais.abinit_input

    def test_runlevel(self, abinit_test_dir):
        abinit_input = load_abinit_input(
            os.path.join(abinit_test_dir, "abinit_inputs"), fname="abinit_input_Si.json"
        )
        ais = AbinitInputSet(abinit_input=abinit_input)
        assert ais.runlevel() == abinit_input.runlevel

    def test_set_structure(self, abinit_test_dir):
        abinit_input = load_abinit_input(
            os.path.join(abinit_test_dir, "abinit_inputs"), fname="abinit_input_Si.json"
        )
        ais = AbinitInputSet(abinit_input=abinit_input)
        struct = abinit_input.structure
        newstruct = struct.copy()
        newstruct.scale_lattice(newstruct.volume * 1.2)
        ais.set_structure(newstruct)
        assert ais.abinit_input.structure == newstruct


@dataclass
class SomeAbinitInputSetGenerator(AbinitInputGenerator):
    calc_type: str = "some_calc"

    param1: int = 1
    param2: Optional[float] = None
    param3: List[int] = field(default_factory=list)
    param4: str = "test_string"

    extra_abivars: dict = field(default_factory=dict)

    restart_from_deps: tuple = (f"{SCF}:WFK|DEN",)

    def get_abinit_input(
        self, structure=None, pseudos=None, prev_outputs=None, **kwargs
    ):
        return AbinitInput(
            structure=structure,
            pseudos=as_pseudo_table(AbinitInputGenerator.pseudos),
        )


class TestAbinitInputSetGenerator:
    def test_from_prev_generator(self):
        saisg1 = SomeAbinitInputSetGenerator(
            param1=2,
            param3=[1, 2, 3],
            extra_abivars={"ecut": 5.0, "nstep": 25, "fake": 1},
        )
        saisg2 = SomeAbinitInputSetGenerator.from_prev_generator(
            prev_input_generator=saisg1,
            param2=1.5,
            param1=10,
            extra_abivars={"ntime": 1, "nstep": None, "fake": 3},
        )
        assert saisg2.calc_type == "some_calc"
        assert saisg2.param1 == 10
        assert saisg2.param2 == 1.5
        assert saisg2.param3 == [1, 2, 3]
        assert saisg2.param4 == "test_string"
        assert saisg2.extra_abivars == {"ecut": 5.0, "fake": 3, "ntime": 1}
        saisg3 = SomeAbinitInputSetGenerator.from_prev_generator(
            prev_input_generator=saisg1, calc_type="new_calc"
        )
        assert saisg3.calc_type == "new_calc"
        with pytest.raises(RuntimeError, match=r"Cannot change pseudos."):
            SomeAbinitInputSetGenerator.from_prev_generator(
                prev_input_generator=saisg1, pseudos="some_pseudo"
            )

    def test_check_format_prev_dirs(self):
        aisg = AbinitInputGenerator()
        prev_outputs = aisg.check_format_prev_dirs(None)
        assert prev_outputs is None
        prev_outputs = aisg.check_format_prev_dirs("/some/path")
        assert prev_outputs == ["/some/path"]
        prev_outputs = aisg.check_format_prev_dirs(Path("/some/path"))
        assert prev_outputs == ["/some/path"]
        prev_outputs = aisg.check_format_prev_dirs(
            ["/some/path", Path("/some/other/path")]
        )
        assert prev_outputs == ["/some/path", "/some/other/path"]

    def test_resolve_dep(self):
        aisg = AbinitInputGenerator()
        with ScratchDir(".") as tmpdir:
            prev_output_dir = os.path.join(tmpdir, "prev_output")
            prev_outdata = os.path.join(prev_output_dir, OUTDIR_NAME)
            makedirs_p(prev_outdata)
            Path(os.path.join(prev_outdata, "out_DEN")).touch()
            irdvars, restart_file = aisg.resolve_dep_exts(
                prev_output_dir, exts=("WFK", "DEN")
            )
            assert irdvars == {"irdden": 1}
            assert restart_file == [(os.path.join(prev_outdata, "out_DEN"), "in_DEN")]
            Path(os.path.join(prev_outdata, "out_WFK")).touch()
            irdvars, restart_file = aisg.resolve_dep_exts(
                prev_output_dir, exts=("WFK", "DEN")
            )
            assert irdvars == {"irdwfk": 1}
            assert restart_file == [(os.path.join(prev_outdata, "out_WFK"), "in_WFK")]
            irdvars, restart_file = aisg.resolve_dep_exts(
                prev_output_dir, exts=("DEN",)
            )
            assert irdvars == {"irdden": 1}
            assert restart_file == [(os.path.join(prev_outdata, "out_DEN"), "in_DEN")]
        with ScratchDir(".") as tmpdir:
            prev_output_dir = os.path.join(tmpdir, "prev_output")
            prev_outdata = os.path.join(prev_output_dir, OUTDIR_NAME)
            makedirs_p(prev_outdata)
            Path(os.path.join(prev_outdata, "out_WFK")).touch()
            irdvars, restart_file = aisg.resolve_dep_exts(
                prev_output_dir, exts=("WFK", "DEN")
            )
            assert irdvars == {"irdwfk": 1}
            assert restart_file == [(os.path.join(prev_outdata, "out_WFK"), "in_WFK")]
            with pytest.raises(
                InitializationError, match=r"Cannot find DDB file to restart from."
            ):
                aisg.resolve_dep_exts(prev_output_dir, exts=("DDB",))
            with pytest.raises(
                InitializationError,
                match=r"Cannot find DDB or DVDB or DKK file to restart from.",
            ):
                aisg.resolve_dep_exts(prev_output_dir, exts=("DDB", "DVDB", "DKK"))
        with ScratchDir(".") as tmpdir:
            prev_output_dir = os.path.join(tmpdir, "prev_output")
            prev_outdata = os.path.join(prev_output_dir, OUTDIR_NAME)
            makedirs_p(prev_outdata)
            Path(os.path.join(prev_outdata, "out_TIM1_DEN")).touch()
            Path(os.path.join(prev_outdata, "out_TIM2_DEN")).touch()
            Path(os.path.join(prev_outdata, "out_TIM15_DEN")).touch()
            irdvars, restart_file = aisg.resolve_dep_exts(
                prev_output_dir, exts=("DEN",)
            )
            assert irdvars == {"irdden": 1}
            assert restart_file == [
                (os.path.join(prev_outdata, "out_TIM15_DEN"), "in_DEN")
            ]
            Path(os.path.join(prev_outdata, "out_DEN")).touch()
            irdvars, restart_file = aisg.resolve_dep_exts(
                prev_output_dir, exts=("DEN",)
            )
            assert irdvars == {"irdden": 1}
            assert restart_file == [(os.path.join(prev_outdata, "out_DEN"), "in_DEN")]

    def test_resolve_deps(self):
        aisg = AbinitInputGenerator()
        with ScratchDir(".") as tmpdir1:
            prev_output_dir1 = os.path.join(tmpdir1, "prev_output")
            prev_outdata1 = os.path.join(prev_output_dir1, OUTDIR_NAME)
            makedirs_p(prev_outdata1)
            den = Path(os.path.join(prev_outdata1, "out_DEN"))
            den.touch()
            with ScratchDir(".") as tmpdir2:
                prev_output_dir2 = os.path.join(tmpdir2, "prev_output")
                prev_outdata2 = os.path.join(prev_output_dir2, OUTDIR_NAME)
                makedirs_p(prev_outdata2)
                wfk = Path(os.path.join(prev_outdata2, "out_WFK"))
                wfk.touch()
                prev_outputs = [prev_output_dir1, prev_output_dir2]
                irdvars, input_files = aisg.resolve_deps(
                    prev_outputs, deps=("any:WFK|DEN",), check_runlevel=False
                )
                assert irdvars == {"irdden": 1, "irdwfk": 1}
                assert (str(wfk), "in_WFK") in input_files
                assert (str(den), "in_DEN") in input_files
                assert len(input_files) == 2

    def test_get_input_set(self, si_structure, mocker):
        with ScratchDir(".") as tmpdir:
            saisg = SomeAbinitInputSetGenerator()
            abinit_input_set = saisg.get_input_set(
                structure=si_structure,
            )
            abinit_input_set.write_input("output1")
            output1 = os.path.join(tmpdir, "output1")

            out_wfk1 = Path(os.path.join(output1, OUTDIR_NAME, "out_WFK"))
            out_wfk1.touch()
            out_den1 = Path(os.path.join(output1, OUTDIR_NAME, "out_DEN"))
            out_den1.touch()

            abinit_input_set.write_input("output1b")
            output1b = os.path.join(tmpdir, "output1b")

            out_den1b = Path(os.path.join(output1b, OUTDIR_NAME, "out_DEN"))
            out_den1b.touch()

            mocker.patch.object(
                atomate2.abinit.sets.base,
                "get_final_structure",
                return_value=si_structure,
            )
            abinit_input_set = saisg.get_input_set(
                structure=si_structure,
                restart_from=output1,
            )
            assert "irdwfk" in abinit_input_set.abinit_input
            assert "irdden" not in abinit_input_set.abinit_input
            assert abinit_input_set.abinit_input["irdwfk"] == 1
            assert abinit_input_set.validate() is True
            abinit_input_set.write_input("output2")
            output2 = os.path.join(tmpdir, "output2")
            in_wfk2 = os.path.join(output2, INDIR_NAME, "in_WFK")
            assert os.path.islink(in_wfk2)
            assert os.readlink(in_wfk2) == str(out_wfk1)

            abinit_input_set = saisg.get_input_set(
                structure=si_structure,
                restart_from=output1b,
            )
            assert "irdwfk" not in abinit_input_set.abinit_input
            assert "irdden" in abinit_input_set.abinit_input
            assert abinit_input_set.abinit_input["irdden"] == 1
            assert abinit_input_set.validate() is True
            abinit_input_set.write_input("output2b")
            output2b = os.path.join(tmpdir, "output2b")
            in_den2b = os.path.join(output2b, INDIR_NAME, "in_DEN")
            assert os.path.islink(in_den2b)
            assert os.readlink(in_den2b) == str(out_den1b)
