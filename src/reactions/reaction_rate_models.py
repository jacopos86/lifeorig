
#
#     Reaction rates factory
#


def build_rate_model(parsed):
    key = parsed.catalyst_or_control.lower()
    # check key
    if key == "2body":
        _, values = parsed.rate_template.split(":")
        A, tfac = map(float, values.split())
        return TwoBodyRate(
            A=A,
            tfac=tfac
        )
    log.error(f"Unknown rate model: {key}")

#
#     Reaction rate models
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

    def __call__(self, temperature):
        return self.A * np.exp(
            self.tfac / temperature
        )