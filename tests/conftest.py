"""Shared pytest configuration: assigns the ``slow``, ``native`` and ``realdata`` markers.

The fast suite is ``pytest -m "not slow" tests``; the full suite (no ``-m``)
remains the release gate. Markers are assigned here, by node id, so that test
files need no edits:

* ``native``: needs a native viewer (PyMOL, VMD or ChimeraX); also ``slow``.
* ``realdata``: needs real structures or the GLUT1 trajectory that are not in
  the repository; also ``slow``.
* ``slow``: took 5 s or more in a single-process run of the full suite
  (October 2026), listed in ``SLOW_TESTS``.
  Re-measure with ``pytest --durations=0`` when tests change.

Every test is xdist-safe: outputs go to ``tmp_path``. Run in parallel with
``pytest -n auto`` (pytest-xdist, in the ``test`` extra).
"""
import pytest

#: Test node ids (parameters stripped) of the slow tests, longest first.
SLOW_TESTS = {
    "tests/test_surface_smoothing.py::test_publication_measurements_identical_with_and_without_smoothing",
    "tests/test_output_inventory.py::test_analyze_single_commands_and_trajectory_write_csv_companions",
    "tests/test_capped_channels.py::test_through_channel_is_unchanged_by_the_option",
    "tests/test_region_definition.py::test_comparison_rejects_unpaired_frames_and_preserves_region_identity",
    "tests/test_cavity_trajectory.py::test_full_trajectory_summary_and_export_are_repeatable",
    "tests/test_summary_statistics_csv.py::test_every_json_interval_statistic_appears_once_in_the_bundle_csvs",
    "tests/test_output_inventory.py::test_publish_default_inventory_is_csv_first_without_html",
    "tests/test_region_radii.py::test_region_compare_refuses_runs_with_different_radius_sets",
    "tests/test_figure_annotations.py::test_publish_data_outputs_identical_with_and_without_annotations",
    "tests/test_region_radii.py::test_region_trajectory_refuses_a_different_set_unless_overridden",
    "tests/test_figure_annotations.py::test_publish_figures_are_clean_by_default_and_annotate_restores_text",
    "tests/test_capped_channels.py::test_exit_leg_casts_are_disjoint_from_the_channel_and_stop_at_bulk",
    "tests/test_radius_sets.py::test_parallel_cavity_trajectory_workers_use_the_chosen_set",
    "tests/test_publish_fallback_and_sticks.py::test_resolved_profile_keeps_extras_and_reports_resolved",
    "tests/test_output_inventory.py::test_network_summary_axis_is_explained_and_counts_residue_contacts",
    "tests/test_auto_enclosure.py::test_publish_records_the_choice_or_the_refusal",
    "tests/test_output_inventory.py::test_chord_colours_follow_the_interaction_names_networks_write",
    "tests/test_auto_enclosure.py::test_static_suite_and_benchmark_default_to_auto",
    "tests/test_region_radii.py::test_prepare_with_an_explicit_set_stores_it_and_trajectory_uses_it",
    "tests/test_region_definition.py::test_region_cli_runs_each_frame_and_retains_region_provenance",
    "tests/test_region_radii.py::test_hydration_cavity_results_uses_or_checks_the_recorded_set",
    "tests/test_publish_fallback_and_sticks.py::test_unresolved_profile_completes_requested_extras_without_a_profile",
    "tests/test_output_policy.py::test_legacy_outputs_reproduce_the_020_file_set_and_numbers",
    "tests/test_output_policy.py::test_exact_inventories_of_default_selected_and_all_modes",
    "tests/test_output_policy.py::test_cli_and_python_give_identical_numbers_and_files",
    "tests/test_output_policy.py::test_renderer_guard_writes_no_scene_or_html_unless_requested",
    "tests/test_output_policy.py::test_default_publish_has_three_viewers_key_figures_and_no_water_files",
    "tests/test_output_policy.py::test_compute_writes_nothing_and_export_matches_the_wrapper",
}

#: Files whose every test needs a native viewer.
NATIVE_FILES = {"test_chimerax_native.py"}
#: Files or classes whose tests need data outside the repository.
REALDATA_FILES = {"test_trajectory_real_gromacs.py", "test_canonical_publish.py"}
REALDATA_CLASSES = {"PDBIntegrationTests"}


def _base_nodeid(item):
    """Node id relative to the repository root, without parameters."""
    nodeid = item.nodeid.split("[", 1)[0]
    return nodeid if nodeid.startswith("tests/") else "tests/" + nodeid.split("tests/", 1)[-1]


def pytest_collection_modifyitems(config, items):
    """Attach the markers described in the module docstring to collected ``items``."""
    for item in items:
        name = item.path.name
        function = getattr(item, "originalname", item.name)
        native = name in NATIVE_FILES or "native" in function
        cls = getattr(item, "cls", None)
        realdata = name in REALDATA_FILES or (cls is not None and cls.__name__ in REALDATA_CLASSES)
        if native:
            item.add_marker(pytest.mark.native)
        if realdata:
            item.add_marker(pytest.mark.realdata)
        if native or realdata or _base_nodeid(item) in SLOW_TESTS:
            item.add_marker(pytest.mark.slow)
