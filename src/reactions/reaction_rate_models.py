

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