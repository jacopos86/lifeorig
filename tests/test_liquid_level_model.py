from src.common.units import Q_
from src.environment.external_drive_params import ExternalDriveParams
from src.environment.liquid_level_model import Liquid, LiquidLevelParams
from src.environment.solvent import SolventData
from src.grids.temporal_grid import TimeGrid


def test_liquid_step_updates_level_and_concentrations():
    liquid = Liquid(
        solvent_data=SolventData(name="test", composition={"CH4": 1.0}),
        liquid_level_params=LiquidLevelParams(base_level=Q_(1.0, "centimeter")),
        external_forces={
            "rainfall": ExternalDriveParams(
                model_type="constant",
                base_level=Q_(0.1, "centimeter / day"),
            ),
            "evaporation": ExternalDriveParams(
                model_type="constant",
                base_level=Q_(0.2, "centimeter / day"),
            ),
        },
        atmospheric_composition={"CO": 1.0},
    )

    liquid.step(dt=Q_(1.0, "day"), time=Q_(0.0, "day"), surface_area=Q_(1.0, "centimeter ** 2"))

    assert liquid.level.to("centimeter").magnitude == 0.9
    assert liquid.concentrations["CH4"] == 0.8 / 0.9
    assert liquid.concentrations["CO"] == 0.1 / 0.9


def test_plot_liquid_level_writes_file_without_changing_level(tmp_path):
    liquid = Liquid(
        solvent_data=SolventData(name="test", composition={"CH4": 1.0}),
        liquid_level_params=LiquidLevelParams(base_level=Q_(1.0, "centimeter")),
        external_forces={
            "rainfall": ExternalDriveParams(
                model_type="constant",
                base_level=Q_(0.1, "centimeter / day"),
            ),
            "evaporation": ExternalDriveParams(
                model_type="constant",
                base_level=Q_(0.0, "centimeter / day"),
            ),
        },
        atmospheric_composition={"CH4": 1.0},
        time_grid=TimeGrid(
            T=Q_(3.0, "day"),
            dt=Q_(1.0, "day"),
            nt=3,
        ),
        working_dir=tmp_path,
    )

    output_file = liquid.plot_liquid_level()

    assert output_file.exists() if hasattr(output_file, "exists") else True
    assert (tmp_path / "liquid_level.pdf").exists()
    assert liquid.level.to("centimeter").magnitude == 1.0
