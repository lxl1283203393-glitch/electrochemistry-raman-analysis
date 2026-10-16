"""Output-contract tests using small, explicitly synthetic spectra."""
from datetime import datetime, timedelta
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd


class FigureOutputTests(unittest.TestCase):
    def test_empty_data_records_all_conditional_skips_without_fake_plots(self):
        from tempfile import TemporaryDirectory
        from optosky_raman.visualize import create_figures

        with TemporaryDirectory() as directory:
            manifest = create_figures([], pd.DataFrame(), pd.DataFrame(), [], {}, Path(directory), {})
            self.assertEqual(len(manifest), 16)
            self.assertTrue(all(row["status"] == "skipped" and row["reason"] for row in manifest))
            self.assertFalse(list(Path(directory).rglob("*.png")))

    def test_applicable_figures_have_three_real_formats_and_unclipped_source_csv(self):
        from tempfile import TemporaryDirectory
        from optosky_raman.visualize import create_figures

        x = np.linspace(200., 3800., 160)
        y = 25 * np.exp(-.5 * ((x - 500) / 50) ** 2) - 3
        records = []
        for index in range(4):
            batch = "B" if index == 3 else "A"
            records.append(SimpleNamespace(
                id=f"{batch}:{index}.csv", path=Path(f"{index}.csv"), batch=batch,
                x=x.copy(), channels={"raw": y + 10, "dark_subtracted": y.copy(), "baseline_subtracted": y / 2},
                stages={"unified": y.copy(), "asls": y.copy(), "arpls": y * .9, "airpls": y * .8},
                metadata={}, state="blank" if index == 0 else "applied_potential",
                potential_v=-1.1 - index * .1, acquired_at=datetime(2026, 1, 1) + timedelta(seconds=index),
                sha256=f"source-{index}", history=[], warnings=[]))
        fits = [{"record_id": records[1].id, "window": [200, 800], "status": "ok", "x": x, "y": y,
                 "fitted": y, "residual": np.zeros_like(y), "baseline": -3 * np.ones_like(y),
                 "components": [{"center_cm1": 500., "area": 100., "fwhm_cm1": 117.7,
                                 "height": 25., "shape": "gaussian", "fitted": y + 3}], "reasons": []}]
        peaks = pd.DataFrame([{"record_id": records[1].id, "position_cm1": 500., "height": 25.}])
        audit = pd.DataFrame([{"record_id": records[2].id, "suspicious": True}])
        differences = pd.DataFrame({"record_id": records[2].id, "previous_id": records[1].id,
                                    "raman_shift_cm1": x, "difference": y * .1})
        with TemporaryDirectory() as directory:
            output = Path(directory)
            manifest = create_figures(records, audit, peaks, fits, {"difference_spectra": differences}, output,
                                      {"experiment": {"reference_electrode": "Ag/AgCl"},
                                       "plotting": {"dpi": 72, "max_fit_panels": 2, "max_annotated_records": 2}})
            self.assertEqual(len(manifest), 16)
            self.assertTrue(all(row["status"] == "generated" for row in manifest), manifest)
            for row in manifest:
                for fmt, filename in row["paths"].items():
                    self.assertEqual(Path(filename).suffix, f".{fmt}")
                    self.assertGreater((output / filename).stat().st_size, 100)
                self.assertEqual(set(row["paths"]), {"png", "pdf", "svg"})
                self.assertGreater((output / row["data_path"]).stat().st_size, 20)
                self.assertIn("Ag/AgCl", row["reference"])
            raw = pd.read_csv(output / manifest[0]["data_path"])
            self.assertLess(raw["intensity"].min(), 0)
            self.assertEqual(set(raw["source_sha256"]), {"source-0", "source-1", "source-2", "source-3"})
            self.assertEqual(len(raw), 4 * len(x))
            self.assertIn("processing", raw.columns)
            self.assertFalse(list(output.glob("*.csv")), "Figure CSVs belong in figures/data, not source root")

    def test_unreliable_components_never_enter_potential_trends(self):
        from tempfile import TemporaryDirectory
        from optosky_raman.visualize import create_figures

        record = SimpleNamespace(id="a", path=Path("a.csv"), batch="A", x=np.array([1., 2., 3.]),
                                 channels={"dark_subtracted": np.array([-1., 0., 1.])}, stages={},
                                 potential_v=-1.2, acquired_at=None, state="applied_potential",
                                 sha256="abc", metadata={}, history=[], warnings=[])
        fits = [{"record_id": "a", "status": "unidentifiable", "components": [
            {"center_cm1": 2., "area": 10., "fwhm_cm1": 1.}]}]
        with TemporaryDirectory() as directory:
            manifest = create_figures([record], pd.DataFrame(), pd.DataFrame(), fits, {}, Path(directory),
                                      {"plotting": {"formats": ["png"], "dpi": 50}})
            trends = [row for row in manifest if row["category"] in {"peak_position_potential", "peak_area_potential", "fwhm_potential"}]
            self.assertEqual(len(trends), 3)
            self.assertTrue(all(row["status"] == "skipped" for row in trends))

    def test_enabled_background_has_before_and_after_traceable_figure_data(self):
        from tempfile import TemporaryDirectory
        from optosky_raman.visualize import create_figures

        def spectrum(name, state, y):
            return SimpleNamespace(id=name, path=Path(name + ".csv"), batch="A", x=np.arange(4.),
                                   channels={"dark_subtracted": np.asarray(y)}, stages={}, potential_v=None,
                                   acquired_at=None, state=state, sha256=name, metadata={}, history=[], warnings=[])
        blank = spectrum("blank", "blank", [2., 2., 2., 2.])
        sample = spectrum("sample", "unknown", [1., 3., 8., 2.])
        sample.stages["background_subtracted"] = np.array([-1., 1., 6., 0.])
        with TemporaryDirectory() as directory:
            output = Path(directory)
            manifest = create_figures([blank, sample], pd.DataFrame(), pd.DataFrame(), [], {}, output,
                                      {"plotting": {"formats": ["png"], "dpi": 50}})
            row = next(item for item in manifest if item["category"] == "blank_comparison")
            frame = pd.read_csv(output / row["data_path"])
            self.assertIn("background_subtracted", set(frame["processing"]))
            self.assertLess(frame.loc[frame.processing == "background_subtracted", "intensity"].min(), 0)

    def test_descending_axis_marker_uses_the_correct_intensity_without_reordering_sources(self):
        from tempfile import TemporaryDirectory
        from optosky_raman.visualize import create_figures

        record = SimpleNamespace(id="a", path=Path("a.csv"), batch="A", x=np.array([3., 2., 1.]),
                                 channels={"dark_subtracted": np.array([-1., 9., 0.])},
                                 stages={"unified": np.array([-1., 9., 0.])}, potential_v=-1.2,
                                 acquired_at=None, state="applied_potential", sha256="a", metadata={}, history=[], warnings=[])
        peaks = pd.DataFrame([{"record_id": "a", "position_cm1": 2., "height": 9.}])
        with TemporaryDirectory() as directory:
            output = Path(directory)
            manifest = create_figures([record], pd.DataFrame(), peaks, [], {}, output,
                                      {"plotting": {"formats": ["png"], "dpi": 50}})
            row = next(item for item in manifest if item["category"] == "candidate_peaks")
            frame = pd.read_csv(output / row["data_path"])
            self.assertIn("plotted_intensity", frame)
            marker = frame.loc[frame.row_type == "candidate_peak"].iloc[0]
            self.assertEqual(marker.plotted_intensity, 9.)
            self.assertEqual(record.x.tolist(), [3., 2., 1.])

    def test_external_background_subtraction_is_not_skipped_without_blank_in_input_records(self):
        from tempfile import TemporaryDirectory
        from optosky_raman.visualize import create_figures

        sample = SimpleNamespace(id="a", path=Path("a.csv"), batch="A", x=np.arange(4.),
                                 channels={"dark_subtracted": np.array([1., 2., 8., 1.])},
                                 stages={"background_subtracted": np.array([-1., 0., 6., -1.])},
                                 potential_v=-1.2, acquired_at=None, state="applied_potential", sha256="a",
                                 metadata={}, warnings=[], history=[{"operation": "background_subtraction",
                                    "scale": 2., "blank_file": "external.csv", "conditions_verified": True}])
        with TemporaryDirectory() as directory:
            output = Path(directory)
            manifest = create_figures([sample], pd.DataFrame(), pd.DataFrame(), [], {}, output,
                                      {"plotting": {"formats": ["png"], "dpi": 50}})
            row = next(item for item in manifest if item["category"] == "blank_comparison")
            self.assertEqual(row["status"], "generated")
            frame = pd.read_csv(output / row["data_path"])
            after = frame.loc[frame.processing == "background_subtracted"]
            self.assertEqual(set(after.blank_file), {"external.csv"})
            self.assertEqual(set(after.background_scale), {2.})


if __name__ == "__main__":
    unittest.main()
