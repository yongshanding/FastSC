
from .bv import get_bv_circuit
from .qaoa import get_qaoa_circuit
from .qgan import get_qgan_circuit
from .ising import get_ising_circuit
from .parallel_cnot import get_parallel_cnot
from .parallel_swap import get_parallel_swap
from .xeb import get_xeb_circuit, get_xeb_iswap_circuit, get_xeb_iswap_barriers_circuit


def get_circuit(numQ, circ_name, dep=0):
    if circ_name=='bv':
        hs = '00101'*(numQ//5 + 1)
        return get_bv_circuit(numQ, hs[:numQ])
    elif circ_name == 'qaoa':
        return get_qaoa_circuit(numQ, 0.5, 1)
    elif circ_name == 'qgan':
        return get_qgan_circuit(numQ)
    elif circ_name == 'ising':
        return get_ising_circuit(numQ)

    elif circ_name == 'parallel_cnot':
        return get_parallel_cnot()
    elif circ_name == 'parallel_swap':
        return get_parallel_swap()
    elif circ_name == 'xeb':
        return get_xeb_circuit(numQ,dep)
    elif circ_name == 'xeb_iswap':
        return get_xeb_iswap_circuit(numQ,dep)
    elif circ_name == 'xeb_iswap_barrier':
        return get_xeb_iswap_barriers_circuit(numQ, dep)
    else:
        print("Circuit name %s not recognized." % circ_name)
        sys.exit()

