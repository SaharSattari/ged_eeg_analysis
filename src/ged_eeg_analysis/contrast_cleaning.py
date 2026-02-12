# REF: Sattari, S., Virji-Babul, N., & Wu, L. C. (2025). 
    #Contrast-based artifact removal enables microstate analysis in ambulatory EEG. 
        #IEEE Transactions on Biomedical Engineering. DOI: [10.1109/TBME.2025.3630112]

import pandas as pd
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


class ArtifactReconstruct:
    """
    Class for performing artifact reconstruction using different cleaning methods.

    Args:
        use_covariance (bool, optional): Whether to use covariance matrix. Defaults to True.
        type (str, optional): The type of cleaning method to use. Can be "ged" for Generalized Eigen Decomposition (Cohen 2022)
        or "cpca" for Contrastive PCA (Abid et al, 2018). Defaults to "ged".
        normalize (bool, optional): Whether to center the data if False or normalize if True. Defaults to False.
        vis_results (bool, optional): Whether to visualize the cleaning results. Defaults to True.
        verbose (bool, optional): Whether to print verbose output. Defaults to False.
    """

    def __init__(
        self,
        use_covariance: bool = True,
        type: str = "ged",
        normalize: bool = False,
        vis_results: bool = True,
        verbose: bool = False,
    ):
        self.use_covariance = use_covariance
        self.type = type
        self.normalize = normalize
        self.verbose = verbose
        self.vis_results = vis_results

    def _check_dimensions(self, data: np.ndarray):
        """
        Check the dimensions of the input data and transpose if necessary.

        Args:
            data (np.ndarray): The input data.

        Returns:
            np.ndarray: The transposed data if necessary.
        """
        return data if data.shape[0] < data.shape[1] else data.T

    def _center_and_normalize(self, data: np.ndarray):
        """
        Center and normalize the input data.

        Args:
            data (np.ndarray): The input data.

        Returns:
            np.ndarray: The centered and normalized data.
        """
        centered = data - np.mean(data, axis=1, keepdims=True)
        return (
            centered / np.std(centered, axis=1, keepdims=True)
            if self.normalize
            else centered
        )

    def _covariance_analysis(self, fg: np.ndarray, bg: np.ndarray, alpha: int):
        """
        Perform covariance analysis for artifact reconstruction.

        Args:
            fg (np.ndarray): The foreground data.
            bg (np.ndarray): The background data.
            alpha (int): The alpha parameter.

        Returns:
            np.ndarray: The eigenvectors obtained from the covariance analysis.
        """
        cov_fg, cov_bg = np.cov(fg, rowvar=True), np.cov(bg, rowvar=True)
        sigma_cov = cov_fg - alpha * cov_bg
        w_cov, v_cov = linalg.eig(sigma_cov)
        sorted_indices_cov = np.argsort(w_cov)[::-1]
        return v_cov[:, sorted_indices_cov]

    def _relative_covariance_analysis(self, fg: np.ndarray, bg: np.ndarray):
        """
        Perform relative covariance analysis for artifact reconstruction.

        Args:
            fg (np.ndarray): The foreground data.
            bg (np.ndarray): The background data.

        Returns:
            Tuple[np.ndarray, np.ndarray]: The eigenvectors and eigenvalues obtained from the relative covariance analysis.
        """
        cov_fg, cov_bg = np.cov(fg, rowvar=True), np.cov(bg, rowvar=True)
        gamma = 0.01
        eigenvalues, _ = np.linalg.eig(cov_bg)
        cov_bg_reg = cov_bg * (1 - gamma) + gamma * np.mean(eigenvalues) * np.eye(
            cov_bg.shape[0]
        )

        try:
            evals, evecs = linalg.eig(cov_fg, cov_bg_reg)
            sidx = np.argsort(evals)[::-1]
            evals = evals[sidx]
            evecs = evecs[:, sidx]

        except np.linalg.LinAlgError:
            print("Matrix is singular and cannot be inverted.")

        return evecs, evals

    def _generate_maps(self, evecs: np.ndarray, fg: np.ndarray):
        """
        Generate artifact maps using the eigenvectors.

        Args:
            evecs (np.ndarray): The eigenvectors.
            fg (np.ndarray): The foreground data.

        Returns:
            np.ndarray: The artifact maps.
        """
        covS = np.cov(fg)
        if evecs.ndim != 2:
            raise ValueError("evecs should be a 2D array")

        Maps = np.zeros_like(evecs)

        for counter in range(evecs.shape[1]):
            Maps[:, counter] = -1 * evecs[:, counter].T @ covS

        return Maps

    def perform_analysis(self, input_fg: np.ndarray, input_bg: np.ndarray, alpha: int):
        """
        Perform artifact reconstruction analysis.

        Args:
            input_fg (np.ndarray): The foreground data.
            input_bg (np.ndarray): The background data.
            alpha (int): The alpha parameter.

        Returns:
            Tuple[Dict[str, np.ndarray], np.ndarray]: The results of the analysis and the artifact maps.
        """
        fg = self._check_dimensions(input_fg)
        bg = self._check_dimensions(input_bg)
        fg_centered = self._center_and_normalize(fg)
        bg_centered = self._center_and_normalize(bg)

        results = {}
        Maps = None

        if self.type == "cpca":
            results["cpca"] = self._covariance_analysis(fg_centered, bg_centered, alpha)
            Maps = self.generate_maps(
                results["cpca"], fg_centered
            )  # TODO should it be fg_centered?

        if self.type == "GED":

            results["ged"] = self._relative_covariance_analysis(
                fg_centered, bg_centered
            )
            evecs_f, _ = results["ged"]

            Maps = self._generate_maps(np.array(evecs_f), fg_centered)

        if self.vis_results == True:

            evecs, evals = results["ged"]
            # plot the eigenspectrum
            plt.figure()
            plt.plot(evals / np.max(evals), "s-", markersize=5, markerfacecolor="k")
            plt.xlim([-0.5, 20.5])
            plt.title("GED eigenvalues")
            plt.xlabel("Component number")
            plt.ylabel("Power ratio (norm-$\lambda$)")

            # component time series
            comp_ts = evecs[:, 0].T @ input_fg
            t = np.arange(10000, 12000)
            plt.figure()
            plt.plot(t, comp_ts[10000:12000])
            plt.title("Time series of first GED component")
            plt.xlabel("Time (samples)")
            plt.ylabel("Amplitude")
            plt.grid(True)
            plt.show()


        #remove first component 

        fg_transformed = evecs.T @ input_fg

        V_pinv = np.linalg.inv(evecs)

        print(V_pinv)

        
        fg_newspace_copy = fg_transformed.copy()

        fg_newspace_copy[:2, :] = 0

        fg_clean = V_pinv.T @ fg_newspace_copy



        return results, Maps, fg_clean

