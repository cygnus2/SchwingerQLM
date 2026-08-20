import numpy as np

from tenpy.algorithms import dmrg
from tenpy.models.model import CouplingMPOModel
from tenpy.models.lattice import Lattice
from tenpy.networks.mps import MPS
from tenpy.networks.site import FermionSite, SpinSite, set_common_charges


class Spin1QLMSchwinger(CouplingMPOModel):
    r"""
    Pedagogical 1D spin-1 quantum-link Schwinger model.

    Unit cell:
        u = 0 : staggered matter fermion at site n
        u = 1 : spin-1 gauge link on bond (n, n+1)

    Hamiltonian:
        H = -t \sum_n ( c_n^\dagger S_n^+ c_{n+1} + h.c. )
            + m \sum_n (-1)^n N_n
            + (g2/2) \sum_n (S_n^z)^2
            + lam_G \sum_n G_n^2
            + lam_b (S_{L-1}^z)^2

    with
        G_n = E_n - E_{n-1} - N_n + bg_n,
        E_n = S_n^z,
        bg_n = (1 - (-1)^n)/2 .
    """

    def init_sites(self, model_params):
        # Matter site: spinless fermion, conserving total fermion number.
        f = FermionSite(conserve='N')

        # Gauge link: spin-1 site, no global Sz conservation.
        s = SpinSite(S=1, conserve=None)

        # Make charge metadata compatible across the mixed unit cell.
        set_common_charges([f, s])

        # Define E^2 = (Sz)^2 as a custom onsite operator on the link.
        s.add_op("E2", s.get_op("Sz Sz"))

        return [f, s]

    def init_lattice(self, model_params):
        L = model_params["L"]  # number of matter sites = number of links in this OBC setup
        bc_MPS = model_params.get("bc_MPS", "finite")
        bc_x = model_params.get("bc_x", "open")
        order = model_params.get("order", "default")

        sites = self.init_sites(model_params)

        # 1D lattice with a 2-site unit cell: [matter, link]
        lat = Lattice(
            Ls=[L],
            unit_cell=sites,
            order=order,
            bc=bc_x,
            bc_MPS=bc_MPS,
            basis=[[2.0]],
            positions=[[0.0], [1.0]],
        )
        return lat

    def init_terms(self, model_params):
        L = model_params["L"]
        t = model_params.get("t", 1.0)
        m = model_params.get("m", 0.0)
        g2 = model_params.get("g2", 1.0)
        lam_G = model_params.get("lam_G", 0.0)
        lam_b = model_params.get("lam_b", 0.0)

        # ------------------------------------------------------------
        # 1) staggered mass term:  m * (-1)^n * N_n
        # matter site is u = 0
        # ------------------------------------------------------------
        mass_array = np.array([m * ((-1) ** n) for n in range(L)], dtype=float)
        self.add_onsite(mass_array, 0, "N")

        # ------------------------------------------------------------
        # 2) electric energy: (g2/2) * E_n^2 = (g2/2) * (Sz)^2
        # link site is u = 1
        # ------------------------------------------------------------
        self.add_onsite(0.5 * g2, 1, "E2")

        # ------------------------------------------------------------
        # 3) gauge-matter hopping:
        #    -t * (c_n^\dagger S_n^+ c_{n+1} + h.c.)
        #
        # In TeNPy add_multi_coupling, each tuple is (opname, dx, u)
        # where dx is the unit-cell displacement.
        # ------------------------------------------------------------
        self.add_multi_coupling(
            -t,
            [
                ("Cd", [0], 0),   # c_n^\dagger     on matter in cell n
                ("Sp", [0], 1),   # S_n^+           on link in cell n
                ("C",  [1], 0),   # c_{n+1}         on matter in cell n+1
            ],
            plus_hc=True,
        )

        # ------------------------------------------------------------
        # 4) Gauss-law penalty: lam_G * sum_n G_n^2
        #
        #    G_n = E_n - E_{n-1} - N_n + bg_n
        #    bg_n = 0 for even n, 1 for odd n
        #
        # We add this as explicit single-site terms using add_onsite_term
        # and explicit fixed couplings using add_local_term.
        # ------------------------------------------------------------
        if lam_G != 0.0:
            for n in range(L):
                bg = 0 if (n % 2 == 0) else 1

                i_m = 2 * n          # MPS index of matter site n
                i_e = 2 * n + 1      # MPS index of link E_n

                # E_n^2 + N_n^2 - 2 E_n N_n + 2 bg E_n - 2 bg N_n
                self.add_onsite_term(lam_G, i_e, "E2")
                self.add_onsite_term(lam_G, i_m, "N")   # since N^2 = N for fermions
                self.add_local_term(-2.0 * lam_G, [("Sz", [n, 1]), ("N", [n, 0])])

                if bg != 0:
                    self.add_onsite_term(2.0 * lam_G * bg, i_e, "Sz")
                    self.add_onsite_term(-2.0 * lam_G * bg, i_m, "N")

                # Terms involving E_{n-1}
                if n > 0:
                    # + E_{n-1}^2
                    i_em1 = 2 * (n - 1) + 1
                    self.add_onsite_term(lam_G, i_em1, "E2")

                    # -2 E_n E_{n-1}
                    self.add_local_term(
                        -2.0 * lam_G,
                        [("Sz", [n - 1, 1]), ("Sz", [n, 1])]
                    )

                    # +2 E_{n-1} N_n
                    self.add_local_term(
                        2.0 * lam_G,
                        [("Sz", [n - 1, 1]), ("N", [n, 0])]
                    )

                    # -2 bg E_{n-1}
                    if bg != 0:
                        self.add_onsite_term(-2.0 * lam_G * bg, i_em1, "Sz")

        # ------------------------------------------------------------
        # 5) right-boundary flux pinning: lam_b * E_{L-1}^2
        # ------------------------------------------------------------
        if lam_b != 0.0:
            i_right_link = 2 * (L - 1) + 1
            self.add_onsite_term(lam_b, i_right_link, "E2")


if __name__ == "__main__":
    L = 12

    model_params = {
        "L": L,
        "bc_MPS": "finite",
        "bc_x": "open",
        "t": 1.0,
        "m": 0.5,
        "g2": 1.0,
        "lam_G": 20.0,
        "lam_b": 20.0,
    }

    model = Spin1QLMSchwinger(model_params)

    # ------------------------------------------------------------
    # Initial product state
    #
    # FermionSite basis is usually:
    #   0 -> empty
    #   1 -> occupied
    #
    # SpinSite(S=1) basis is usually:
    #   0 -> m=-1
    #   1 -> m= 0
    #   2 -> m=+1
    #
    # So choose staggered bare vacuum and zero electric field.
    # ------------------------------------------------------------
    product_state = []
    for n in range(L):
        product_state.append(0 if (n % 2 == 0) else 1)  # matter site
        product_state.append(1)                         # link site, E = 0

    psi = MPS.from_product_state(model.lat.mps_sites(), product_state, bc="finite", unit_cell_width=2)

    dmrg_params = {
        "mixer": True,
        "trunc_params": {
            "chi_max": 100,
            "svd_min": 1e-10,
        },
        "max_E_err": 1e-8,
        "max_sweeps": 20,
        "combine": True
    }

    info = dmrg.run(psi, model, dmrg_params)
    #print(info)
    ss = info["sweep_statistics"]
    print("Final energy         =", info["E"])
    print("Number of sweeps     =", ss["sweep"][-1])
    print("Total updates        =", ss["N_updates"][-1])
    print("Total time (s)       =", ss["time"][-1])
    print("Peak max_chi         =", max(ss["max_chi"]))
    print("Final sweep max_chi  =", ss["max_chi"][-1])
    print("Final max_trunc_err  =", ss["max_trunc_err"][-1])
    print("Final max_E_trunc    =", ss["max_E_trunc"][-1])
    print("Final norm_err       =", ss["norm_err"][-1])
    
    print("Final MPS chi list   =", psi.chi)
    print("Final MPS max chi    =", max(psi.chi))
