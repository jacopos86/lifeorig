from dataclasses import dataclass
from pathlib import Path
from src.global_parameters.global_variables import VPL_ARCHEAN_REACTIONS
from src.parallelization.array_backend import xp
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
            reaction_id=parsed.reaction_id,
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
            return TwoTermArrheniusRate(
                A1=6.34e-12,
                tfac1=-4000.0,
                A2=1.46e-9,
                tfac2=-9650.0
            )
        
        #   reaction : VPL0016
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0016"
        ):
            return TwoTermDensityRate(
                A1=3.0e-13,
                tfac1=460.0,
                A2=2.1e-33,
                tfac2=920.0
            )

        #   reaction : VPL0018
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0018"
        ):
            return DensityPowerLawRate(
                A=6.0e-34,
                T_ref=300.0,
                temperature_power=-2.6
            )

        #   reaction : VPL0032
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0032"
        ):
            return COOxidationFallOffRate(
                k0=1.5e-13,
                Ainf=2.1e9,
                T_ref=300.0,
                temperature_power=6.1,
                broad_coeff=0.6
            )

        #   reaction : VPL0033
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0033"
        ):
            return DensityScaledArrheniusRate(
                A=1.6e-32,
                tfac=-2184.0
            )

        #   reaction : VPL0034
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0034"
        ):
            return DensityScaledArrheniusRate(
                A=1.0e-33,
                tfac=-629.0
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
                tfac0=-500.0,
                planet_spec_factors={
                    "MARS": 2.5
                }
            )

        #   reaction: VPL0074
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0074"
        ):
            return ModifiedArrheniusRate(
                alpha=4.9e-15,
                beta=4.4,
                gamma=2450.0,
                T_ref=298.0
            )

        #   reaction: VPL0075
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0075"
        ):
            return ActivatedThreeBodyFalloffRate(
                A0=1.2e-31,
                Ainf=2.4e-10,
                n0=1.17,
                ninf=0.41,
                tfac0=-210.0
            )

        #   reaction: VPL0076
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0076"
        ):
            return DensityScaledArrheniusRate(
                A=1.25e-32,
                tfac=0.0,
                density_power=1.0
            )

        #   reaction: VPL0090
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0090"
        ):
            return TwoBodyRate(
                A=3.0e-11,
                tfac=-500.0
            )

        #   reaction: VPL0101
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0101"
        ):
            return ModifiedArrheniusRate(
                alpha=8.54e-12,
                beta=1.5,
                gamma=2920.0,
                T_ref=300.0
            )

        #   reaction: VPL0107
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0107"
        ):
            return TwoBodyRate(
                A=0.0,
                tfac=0.0
            )

        #   reaction: VPL0108
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0108"
        ):
            return HNO3OHRate()

        #   reaction: VPL0118
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0118"
        ):
            return DensityScaledArrheniusRate(
                A=6.0e-31,
                tfac=0.0,
                density_power=1.0
            )

        #   reaction: VPL0142
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0142"
        ):
            return DensityScaledArrheniusRate(
                A=1.87e-33,
                tfac=-206.0,
                density_power=1.0
            )

        #   reaction: VPL0146
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0146"
        ):
            return CappedDensityRate(
                A=2.5e-30 / 5.0e1,
                maximum_rate=5.0e-11
            )

        #   reaction: VPL0147
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0147"
        ):
            return CappedDensityRate(
                A=2.5e-30 / 1.0e1,
                maximum_rate=5.0e-11
            )

        # reaction: VPL0148
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0148"
        ):
            return CappedDensityRate(
                A=2.5e-30 / 5.0e1,
                maximum_rate=5.0e-11
            )

        # reaction: VPL0149
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0149"
        ):
            return CappedDensityRate(
                A=2.5e-30 / 1.0e1,
                maximum_rate=5.0e-11
            )

        # reaction: VPL0190
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0190"
        ):
            return DensityScaledArrheniusRate(
                A=4.0e-33,
                tfac=-1780.0,
                density_power=1.0
            )

        # reaction: VPL0193
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0193"
        ):
            return DensityScaledArrheniusRate(
                A=8.3e-33,
                tfac=0.0,
                density_power=1.0
            )

        # reaction: VPL0209
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0209"
        ):
            return PressureScaledRate(
                A=1.5e-13,
                pressure_factor=0.6
            )

        # reaction: VPL0212
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0212"
        ):
            return ModifiedArrheniusRate(
                alpha=1.9e-14,
                beta=3.97,
                gamma=580.0,
                T_ref=300.0
            )

        # reaction: VPL0217
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0217"
        ):
            return DensityScaledArrheniusRate(
                A=2.0e-33,
                tfac=-850.0,
                density_power=1.0
            )

        # reaction: VPL0222
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0222"
        ):
            return TemperatureDependentLindemannRate(
                A0=5.5e-23,
                temperature_power0=-2.0,
                tfac0=-1040.0,
                Ainf=1.5e-13,
                temperature_power_inf=0.0,
                tfac_inf=-440.0
            )

        # reaction: VPL0239
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0239"
        ):
            return TemperatureDependentLindemannRate(
                A0=8.75e-31,
                temperature_power0=0.0,
                tfac0=524.0,
                Ainf=8.3e-11,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0246
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0246"
        ):
            return TemperatureDependentLindemannRate(
                A0=8.75e-31,
                temperature_power0=0.0,
                tfac0=524.0,
                Ainf=8.3e-11,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0256
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0256"
        ):
            return TemperatureDependentLindemannRate(
                A0=3.1e-30,
                temperature_power0=0.0,
                tfac0=457.0,
                Ainf=1.5e-10,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0267
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0267"
        ):
            return ModifiedArrheniusRate(
                alpha=2.3e-12,
                beta=-0.7,
                gamma=0.0,
                T_ref=300.0
            )

        # reaction: VPL0272
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0272"
        ):
            return TemperatureDependentLindemannRate(
                A0=1.0e-28,
                temperature_power0=0.0,
                tfac0=0.0,
                Ainf=1.0e-15,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0275
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0275"
        ):
            return DensityScaledArrheniusRate(
                A=1.4e-32,
                tfac=-3000.0,
                density_power=1.0
            )

        # reaction: VPL0288
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0288"
        ):
            return TemperatureDependentLindemannRate(
                A0=1.26e-18,
                temperature_power0=-3.1,
                tfac0=-721.0,
                Ainf=3.0e-10,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0299
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0299"
        ):
            return TemperatureDependentLindemannRate(
                A0=2.6e-31,
                temperature_power0=0.0,
                tfac0=0.0,
                Ainf=3.8e-11,
                temperature_power_inf=0.0,
                tfac_inf=-1374.0
            )

        # reaction: VPL0302
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0302"
        ):
            return ModifiedArrheniusRate(
                alpha=2.4e-24,
                beta=4.02,
                gamma=2754.0,
                T_ref=1.0
            )

        # reaction: VPL0306
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0306"
        ):
            return TemperatureDependentLindemannRate(
                A0=2.15e-29,
                temperature_power0=0.0,
                tfac0=-349.0,
                Ainf=4.95e-11,
                temperature_power_inf=0.0,
                tfac_inf=-1051.0
            )

        # reaction: VPL0308
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0308"
        ):
            return CappedArrheniusRate(
                A=2.5e-11,
                tfac=200.0,
                maximum_rate=1.7e-10
            )

        # reaction: VPL0316
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0316"
        ):
            return ModifiedArrheniusRate(
                alpha=3.25e-11,
                beta=-0.5,
                gamma=0.0,
                T_ref=1.0
            )

        # reaction: VPL0320
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0320"
        ):
            return TemperatureDependentLindemannRate(
                A0=5.8e-31,
                temperature_power0=0.0,
                tfac0=1258.0,
                Ainf=1.4e-12,
                temperature_power_inf=0.0,
                tfac_inf=388.0
            )

        # reaction: VPL0336
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0336"
        ):
            return TemperatureDependentLindemannRate(
                A0=2.519e-16,
                temperature_power0=-2.458,
                tfac0=0.0,
                Ainf=8.12e-10,
                temperature_power_inf=-0.5,
                tfac_inf=0.0
            )

        # reaction: VPL0338
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0338"
        ):
            return TwoTermArrheniusRate(
                A1=1.6e-11,
                tfac1=-2900.0,
                A2=2.2e-11,
                tfac2=-2200.0
            )

        # reaction: VPL0351
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0351"
        ):
            return TemperatureDependentLindemannRate(
                A0=1.3e-22,
                temperature_power0=0.0,
                tfac0=0.0,
                Ainf=1.2e-10,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0354
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0354"
        ):
            return TemperatureDependentLindemannRate(
                A0=2.15e-29,
                temperature_power0=0.0,
                tfac0=-349.0,
                Ainf=4.95e-11,
                temperature_power_inf=0.0,
                tfac_inf=-1051.0
            )

        # reaction: VPL0358
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0358"
        ):
            return CappedArrheniusRate(
                A=5.5e-11,
                tfac=173.0,
                maximum_rate=3.55e-10
            )

        # reaction: VPL0359
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0359"
        ):
            return CappedArrheniusRate(
                A=5.5e-11,
                tfac=173.0,
                maximum_rate=3.55e-10
            )

        # reaction: VPL0360
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0360"
        ):
            return TemperatureDependentLindemannRate(
                A0=8.0e-24,
                temperature_power0=-2.0,
                tfac0=-1225.0,
                Ainf=9.7e-13,
                temperature_power_inf=0.0,
                tfac_inf=-1550.0
            )

        # reaction: VPL0361
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0361"
        ):
            return TemperatureDependentLindemannRate(
                A0=8.0e-24,
                temperature_power0=-2.0,
                tfac0=-1225.0,
                Ainf=1.4e-11,
                temperature_power_inf=0.0,
                tfac_inf=-1000.0
            )

        # reaction: VPL0363
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0363"
        ):
            return PressureSuppressedLindemannRate(
                k0=1.9e-27,
                k_inf=2.5e-11
            )

        # reaction: VPL0364
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0364"
        ):
            return TemperatureDependentLindemannRate(
                A0=1.0e-28,
                temperature_power0=0.0,
                tfac0=0.0,
                Ainf=1.0e-11,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0368
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0368"
        ):
            return CappedArrheniusRate(
                A=1.75e-10,
                tfac=61.0,
                maximum_rate=5.3e-10
            )

        # reaction: VPL0369
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0369"
        ):
            return TemperatureDependentLindemannRate(
                A0=3.8e-25,
                temperature_power0=0.0,
                tfac0=0.0,
                Ainf=2.2e-12,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0370
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0370"
        ):
            return TemperatureDependentLindemannRate(
                A0=8.0e-24,
                temperature_power0=-2.0,
                tfac0=-1225.0,
                Ainf=9.7e-12,
                temperature_power_inf=0.0,
                tfac_inf=-1550.0
            )

        # reaction: VPL0371
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0371"
        ):
            return TemperatureDependentLindemannRate(
                A0=8.0e-24,
                temperature_power0=-2.0,
                tfac0=-1225.0,
                Ainf=9.7e-12,
                temperature_power_inf=0.0,
                tfac_inf=-1550.0
            )

        # reaction: VPL0372
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0372"
        ):
            return TemperatureDependentLindemannRate(
                A0=1.7e-26,
                temperature_power0=0.0,
                tfac0=0.0,
                Ainf=1.5e-10,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0373
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0373"
        ):
            return TemperatureDependentLindemannRate(
                A0=1.7e-26,
                temperature_power0=0.0,
                tfac0=0.0,
                Ainf=1.5e-10,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0374
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0374"
        ):
            return TemperatureDependentLindemannRate(
                A0=1.7e-26,
                temperature_power0=0.0,
                tfac0=0.0,
                Ainf=1.5e-10,
                temperature_power_inf=0.0,
                tfac_inf=0.0
            )

        # reaction: VPL0378
        #
        if weird_key == (
            VPL_ARCHEAN_REACTIONS,
            "VPL0378"
        ):
            return TemperatureDependentLindemannRate(
                A0=3.8e-25,
                temperature_power0=0.0,
                tfac0=0.0,
                Ainf=3.7e-12,
                temperature_power_inf=0.0,
                tfac_inf=0.0
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
    planet_spec_factors: dict[str, float] | None = None
    def __call__(self, context, *, collider_density):
        T = context.temperature.to("kelvin").magnitude
        D_M = collider_density.to("1 / centimeter**3").magnitude
        # planet specific factor
        planet_factor = 1.0
        if self.planet_spec_factors and context.planet:
            planet_factor = self.planet_spec_factors.get(
                context.planet.upper(),
                1.0
            )
        k0 = (
            self.A0 * planet_factor
            * xp.exp(self.tfac0 / T)
            * (300.0 / T) ** self.n0
        )
        kinf = (
            self.Ainf * (300.0 / T) ** self.ninf
        )
        red_pressure = (
            k0 * D_M / kinf
        )
        broad_exp = (
            1.0 / (1.0 + xp.log10(red_pressure) ** 2)
        )
        return (
            k0 * D_M / (1.0 + red_pressure)
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
class TemperatureDependentLindemannRate:
    A0: float
    temperature_power0: float
    tfac0: float
    Ainf: float
    temperature_power_inf: float
    tfac_inf: float
    def __call__(self, context, *, collider_density):
        T = context.temperature.to("kelvin").magnitude
        D_M = collider_density.to("1 / centimeter**3").magnitude
        k0 = (
            self.A0 * T ** self.temperature_power0
            * xp.exp(self.tfac0 / T)
        )
        k_inf = (
            self.Ainf
            * T ** self.temperature_power_inf
            * xp.exp(self.tfac_inf / T)
        )
        return (
            k0 * k_inf * D_M
            / (k0 * D_M + k_inf)
        )

@dataclass
class PressureSuppressedLindemannRate:
    k0: float
    k_inf: float
    def __call__(self, context, *, collider_density):
        D_M = collider_density.to("1 / centimeter**3").magnitude
        return (
            self.k_inf ** 2
            / (self.k0 * D_M + self.k_inf)
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
class TwoTermArrheniusRate:
    A1: float
    tfac1: float
    A2: float
    tfac2: float
    def __call__(self, context):
        T = context.temperature.to("kelvin").magnitude
        return (
            self.A1 * xp.exp(self.tfac1 / T)
            + self.A2 * xp.exp(self.tfac2 / T)
        )

@dataclass
class CappedArrheniusRate:
    A: float
    tfac: float
    maximum_rate: float
    def __call__(self, context):
        T = context.temperature.to("kelvin").magnitude
        return xp.minimum(
            self.A * xp.exp(self.tfac / T),
            self.maximum_rate
        )

@dataclass
class DensityScaledArrheniusRate:
    A: float
    tfac: float
    density_power: float = 1.0
    def __call__(self, context, *, collider_density):
        T = context.temperature.to("kelvin").magnitude
        D_M = collider_density.to("1 / centimeter**3").magnitude
        return (
            self.A * xp.exp(self.tfac / T)
            * D_M ** self.density_power
        )

@dataclass
class DensityPowerLawRate:
    A: float
    T_ref: float
    temperature_power: float
    def __call__(self, context, *, collider_density):
        T = context.temperature.to("kelvin").magnitude
        D_M = collider_density.to("1 / centimeter**3").magnitude
        return (
            self.A 
            * (T / self.T_ref) ** self.temperature_power
            * D_M
        )

@dataclass
class TwoTermDensityRate:
    A1: float
    tfac1: float
    A2: float
    tfac2: float
    def __call__(self, context, *, collider_density):
        T = context.temperature.to("kelvin").magnitude
        D_M = collider_density.to("1 / centimeter**3").magnitude
        return (
            self.A1 * xp.exp(self.tfac1 / T)
            + self.A2 * xp.exp(self.tfac2 / T) * D_M
        )

@dataclass
class PressureScaledRate:
    A: float
    pressure_factor: float
    def __call__(self, context):
        P_ATM = context.pressure.to("atmosphere").magnitude
        return self.A * (
            1.0 +
            self.pressure_factor * P_ATM
        )

@dataclass
class COOxidationFallOffRate:
    k0: float
    Ainf: float
    T_ref: float
    temperature_power: float
    broad_coeff: float
    def __call__(self, context, *, collider_density):
        T = context.temperature.to("kelvin").magnitude
        D_M = collider_density.to("1 / centimeter**3").magnitude
        k_inf = (
            self.Ainf
            * (T / self.T_ref) ** self.temperature_power
        )
        red_pressure = self.k0 * D_M / k_inf
        broad_exp = 1.0 / (
            1.0 + xp.log10(red_pressure) ** 2
        )
        return (
            self.k0
            / (1.0 + red_pressure)
            * self.broad_coeff ** broad_exp
        )

@dataclass
class HNO3OHRate:
    def __call__(self, context, *, collider_density):
        T = context.temperature.to("kelvin").magnitude
        D_M = collider_density.to("1 / centimeter**3").magnitude
        k0 = 2.4e-14 * xp.exp(460.0 / T)
        k2 = 2.7e-17 * xp.exp(2199.0 / T)
        k3_M = (
            6.5e-34
            * xp.exp(1335.0 / T) * D_M
        )
        return k0 + k3_M / (1.0 + k3_M / k2)

@dataclass
class CappedDensityRate:
    A: float
    maximum_rate: float
    def __call__(self, context, *, collider_density):
        D_M = collider_density.to("1 / centimeter**3").magnitude
        return xp.minimum(
            self.maximum_rate,
            self.A * D_M
        )

@dataclass
class PiecewiseModifiedArrheniusRate:
    branches: tuple

@dataclass
class PhotolysisRate:
    reaction_id: str
    photo_type: str
    source_data: str
    def __call__(self, context):
        if context.photolysis_rates is None:
            log.error(
                "Photolysis rates are missing from the reaction context"
            )
        if self.reaction_id not in context.photolysis_rates:
            raise KeyError(
                f"Missing photolysis rate for {self.reaction_id}"
            )
        return context.photolysis_rates[self.reaction_id]