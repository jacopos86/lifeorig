from dataclasses import dataclass
from multiprocessing import context
from pathlib import Path
from src.global_parameters.global_variables import VPL_ARCHEAN_REACTIONS
from src.utilities.logging_module import log

#
#     Reaction rates factory
#

def build_rate_model(parsed):
    key = parsed.catalyst_or_control.lower()
    # check key
    if key == "2body":
        _, values = parsed.rate_template.split(":", 1)
        values = values.replace("D", "E").replace("d", "e")
        parameters = list(map(float, values.split()))
        if not 1 <= len(parameters) <= 3:
            log.error(
                f"Invalid 2BODY parameters for {parsed.reaction_id}: "
                f"{parsed.rate_template}"
            )
        A = parameters[0]
        tfac = parameters[1] if len(parameters) >= 2 else 0.0
        return TwoBodyRate(
            A=A,
            tfac=tfac
        )
    if key == "3body":
        _, values = parsed.rate_template.split(":", 1)
        values = values.replace("D", "E").replace("d", "e")
        A0, Ainf, n0, ninf = map(float, values.split())
        return ThreeBodyFalloffRate(
            A0=A0,
            Ainf=Ainf,
            n0=n0,
            ninf=ninf
        )
    if key == "three_body":
        model_name, parameter_text = parsed.rate_template.split("(", 1)
        if model_name.strip().lower() != "lindemann":
            log.error(
                f"Unknown three-body model for {parsed.reaction_id}: "
                f"{model_name}"
            )
        parameter_text = parameter_text.rstrip(")")
        parameters = dict(
            item.strip().split("=", 1)
            for item in parameter_text.split(",")
        )
        k_inf = float(parameters["k_inf"])
        k0 = float(parameters["k0"])
        T_ref = float(parameters["T_ref"].rstrip("Kk"))
        return LindemannRate(
            k0=k0,
            k_inf=k_inf,
            T_ref=T_ref
        )
    if key == "thermal":
        model_name, parameter_text = parsed.rate_template.split("(", 1)
        if model_name.strip().lower() == "piecewise_modified_arrhenius":
            parameter_text = parameter_text.rstrip(")")
            branches = []
            for branch_text in parameter_text.split(";"):
                parameters = dict(
                    item.strip().split("=", 1)
                    for item in branch_text.split(",")
                )
                branches.append(
                    ModifiedArrheniusRate(
                        alpha=float(parameters["alpha"]),
                        beta=float(parameters["beta"]),
                        gamma=float(parameters["gamma"]),
                        T_ref=300.0
                    )
                )
            return PiecewiseModifiedArrheniusRate(
                branches=tuple(branches)
            )
        if model_name.strip().lower() != "modified_arrhenius":
            log.error(
                f"Unknown thermal model for {parsed.reaction_id}: "
                f"{model_name}"
            )
        parameter_text = parameter_text.rstrip(")")
        parameters = dict(
            item.strip().split("=", 1)
            for item in parameter_text.split(",")
        )
        return ModifiedArrheniusRate(
            alpha=float(parameters["alpha"]),
            beta=float(parameters["beta"]),
            gamma=float(parameters["gamma"]),
            T_ref=float(parameters["T_ref"].rstrip("Kk"))
        )
    if key == "photo":
        photo_type, source_data = parsed.rate_template.split(":", 1)
        return PhotolysisRate(
            photo_type=photo_type.strip(),
            source_data=source_data.strip()
        )
    if key == "weird":
        if parsed.source_file is None:
            log.error(
                f"Missing source file for reaction "
                f"{parsed.reaction_id}"
            )
        weird_key = (
            Path(parsed.source_file).resolve(),
            parsed.reaction_id,
        )
        
        #   reaction : VPL0004
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0004"
        ):
            return ModifiedArrheniusRate(
                alpha=1.34e-15,
                beta=6.52,
                gamma=1460.0,
                T_ref=298.0
            )
        
        #   reaction (: VPL0016
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0016"
        ):
            return TwoTermDensityRate(
                A1=2.3e-13,
                tfac1=590.0,
                A2=1.7e-33,
                tfac2=1000.0
            )

        #   reaction : VPL0018
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0018"
        ):
            return DensityScaledArrheniusRate(
                A=9.46e-34,
                tfac=480.0
            )

        #   reaction : VPL0032
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0032"
        ):
            return PressureScaledRate(
                A=1.5e-13,
                pressure_factor=0.6
            )

        #   reaction : VPL0033
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0033"
        ):
            return DensityScaledArrheniusRate(
                A=2.2e-33,
                tfac=-1780.0
            )

        #   reaction : VPL0034
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0034"
        ):
            return DensityScaledArrheniusRate(
                A=1.4e-34,
                tfac=-100.0
            )

        #   reaction : VPL0043
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0043"
        ):
            return ModifiedArrheniusRate(
                alpha=2.14e-12,
                beta=1.62,
                gamma=1090.0,
                T_ref=298.0
            )

        #   reaction : VPL0044
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0044"
        ):
            return DensityPowerLawRate(
                A=8.85e-33,
                T_ref=287.0,
                temperature_power=-0.6
            )

        #   reaction : VPL0047
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0047"
        ):
            return DensityPowerLawRate(
                A=6.9e-31,
                T_ref=298.0,
                temperature_power=-2.0
            )

        #   reaction: VPL0069
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0069"
        ):
            return ActivatedThreeBodyFalloffRate(
                A0=1.17e-25,
                Ainf=3.0e-11,
                n0=3.75,
                ninf=1.0,
                tfac0=-500.0
            )
        
        raise NotImplementedError(
            f"WEIRD rate not implemented for "
            f"source={parsed.source_file}, "
            f"reaction_id={parsed.reaction_id}, "
            f"template={parsed.rate_template}"
        )
    
    log.error(f"Unknown rate model: {key}")

#
#     Reaction rate models -> generic
#

def reaction_rate_gas_phase_model(Delta, x0, xgr, sig, E_a, T):
    ''' Define reaction rates valid 
    for gas phase reactions'''
    N = len(xgr)
    rr = np.zeros(N)
    for i in range(N):
        rr[i] = Delta * exp(-(xgr[i] - x0)**2 / (2*sig**2)) * exp(-E_a / (R*T))
    return rr

def reaction_rate_surface_catalyst_model():
    ''' Define reaction rates with surface catalysts '''
    pass

#
#     Reaction rates -> specific to networks
#

@dataclass
class TwoBodyRate:
    A: float
    tfac: float
    def __call__(self, context):
        T = context.temperature.to("kelvin").magnitude
        return self.A * xp.exp(
            self.tfac / T
        )

@dataclass
class ThreeBodyFalloffRate:
    A0: float
    Ainf: float
    n0: float
    ninf: float
    def __call__(self, context, *, collider_density):
        T = context.temperature.to("kelvin").magnitude
        D_M = collider_density.to("1 / centimeter**3").magnitude
        k0 = self.A0 * (300.0 / T) ** self.n0
        kinf = self.Ainf * (300.0 / T) ** self.ninf
        red_pressure = k0 * D_M / kinf
        broad_exponent = 1.0 / (
            1.0 + xp.log10(red_pressure) ** 2
        )
        return (
            k0 * D_M / (1.0 + red_pressure) * 0.6 ** broad_exponent
        )

@dataclass
class ActivatedThreeBodyFalloffRate:
    A0: float
    Ainf: float
    n0: float
    ninf: float
    tfac0: float
    broad_coeff: float = 0.6
    def __call__(self, context, *, collider_density):
        T = context.temperature.to("kelvin").magnitude
        D_M = collider_density.to("1 / centimeter**3").magnitude
        if self.collision_efficiencies:
            for species_name, efficiency in (
                self.collision_efficiencies.items()
            ):
                species_index = context.species_to_index.get(
                    species_name
                )
                # if species is absent does not
                #  contribute
                if species_index is None:
                    continue
                species_density = (
                    species_D[
                        ...,
                        species_index
                    ]
                )
                # The baseline total density already includes this
                # species once, so add only its excess efficiency.
                effective_density = (
                    effective_density
                    + (efficiency - 1.0) * species_density
                )
        k0 = (
            self.A0 * xp.exp(self.tfac0 / T)
            * (300.0 / T) ** self.n0
        )
        kinf = (
            self.Ainf * (300.0 / T) ** self.ninf
        )
        red_pressure = (
            k0 * effective_density / kinf
        )
        broad_exp = (
            1.0 / (1.0 + xp.log10(red_pressure))
        )
        return (
            k0 * effective_density / (1.0 + red_pressure)
            * self.broad_coeff ** broad_exp
        )

@dataclass
class LindemannRate:
    k0: float
    k_inf: float
    T_ref: float
    def __call__(self, context):
        D = context.number_density.to("1 / centimeter**3").magnitude
        reduced_pressure = (
            self.k0 * D / self.k_inf
        )
        return (
            self.k0 * D / (1.0 + reduced_pressure)
        )

@dataclass
class ModifiedArrheniusRate:
    alpha: float
    beta: float
    gamma: float
    T_ref: float
    def __call__(self, context):
        T = context.temperature.to("kelvin").magnitude
        return (
            self.alpha * (T / self.T_ref) ** self.beta
            * xp.exp(-self.gamma / T)
        )

@dataclass
class DensityScaledArrheniusRate:
    A: float
    tfac: float
    density_power: float = 1.0
    def __call__(self, context):
        T = context.temperature
        D = context.number_density
        return (
            self.A * xp.exp(self.tfac / T)
            * D ** self.density_power
        )

@dataclass
class DensityPowerLawRate:
    A: float
    T_ref: float
    temperature_power: float
    def __call__(self, context):
        return (
            self.A * (
                context.temperature / self.T_ref
            ) ** self.temperature_power
            * context.number_density
        )

@dataclass
class TwoTermDensityRate:
    A1: float
    tfac1: float
    A2: float
    tfac2: float
    def __call__(self, context):
        T = context.temperature
        D = context.number_density
        return (
            self.A1 * xp.exp(self.tfac1 / T)
            + self.A2 * xp.exp(self.tfac2 / T) * D
        )

@dataclass
class PressureScaledRate:
    A: float
    pressure_factor: float
    def __call__(self, context):
        return self.A * (
            1.0 +
            self.pressure_factor * context.pressure_atm
        )

@dataclass
class PiecewiseModifiedArrheniusRate:
    branches: tuple

@dataclass
class PhotolysisRate:
    photo_type: str
    source_data: str