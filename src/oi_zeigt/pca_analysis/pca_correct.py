from __future__ import print_function

import copy
import gc
import math
import os
import pickle
import subprocess
import time
import uuid

import numpy as np
import numpy.ma as ma
import pandas as pd
import pgutils
import pyclass
from sicparse import OptionParser

# import matplotlib.pyplot as plt
from sklearn import preprocessing

from .pca_errors import (
    InvalidOption,
    MissingMandatoryOption,
    MultipleScienceSourcesFound,
    NoScienceSourceFound,
)
from .pca_utilities import (
    add_common_options,
    add_correct_options,
    apply_pca_exclude_range,
    create_index,
    create_pickle_file_name,
    derive_channel_from_velocity,
    find_lines,
    get_and_check_science_sources,
    get_exclude_channels,
    get_line_windows,
    load_config_file,
    prepare_spectrum,
    refill_pca_exclude_range,
    set_line_windows,
)


def plot_pca_decomposition(
    fig=None,
    pca=None,
    decomposition=None,
    number_components=5,
    x_axis=None,
    x_axis_unit="km/s",
    original_spectra=None,
    corrected_spectra=None,
    cutoff=None,
    decomposition_type="pca",
    channels_exclude=None,
    bad_channels_orig=None,
    good_channels_orig=None,
    bad_channels=None,
    good_channels=None,
    norms=None,
    cut_coefficients=0.1,
    global_structure=None,
    plot_width=0.15,
    padding=0.02,
):
    # Plot the results
    import matplotlib.ticker as ticker

    extra_components = 0  # counts spectra shown in middle spectrum
    if original_spectra is not None:
        original_mean = ma.mean(original_spectra, 0)
        if type(original_mean) == np.float64:
            original_mean = original_spectra
        extra_components += 1
    if corrected_spectra is not None:
        corrected_mean = ma.mean(corrected_spectra, 0)
        if type(corrected_mean) == np.float64:
            corrected_mean = corrected_spectra
        extra_components += 1
    spectral_height = 0.85 / (number_components + extra_components)
    i_drawn = 0
    min_value = np.min([np.min(original_mean), np.min(corrected_mean)])
    min_value = min_value - 0.1 * min_value
    max_value = np.max([np.max(corrected_mean), np.max(corrected_mean)])
    max_value = max_value + 0.1 * max_value
    if original_spectra is not None:
        ax = fig.add_axes(
            [padding, 0.90 - spectral_height, plot_width, spectral_height]
        )
        # ax.yaxis.set_major_formatter(ticker.NullFormatter())
        if np.any(bad_channels_orig):
            try:
                original_mean[bad_channels_orig] = 0
            except IndexError:
                pass
        # ax.xaxis.set_major_locator(ticker.MultipleLocator(len(x_axis)))
#        for i in range(len(x_axis)):
#            print(f"{x_axis[i]},  {original_mean[i]}")
        ax.plot(x_axis, original_mean, lw=1, color="red")
        title = "original mean"
        if np.any(bad_channels_orig):
            title += " Blanked values present"
        ax.text(0.03, 0.94, title, transform=ax.transAxes, ha="left", va="top")
        ax.set_ylim(min_value, max_value)
        for line_ in ax.get_xticklines() + ax.get_yticklines():
            line_.set_markersize(2)
        if i_drawn == 0:
            ax.set_title("PCA decomposition into Eigenspectra")
        i_drawn += 1
        ax.axhline(0, color="gray")

    if corrected_spectra is not None:
        ax = fig.add_axes(
            [padding, 0.90 - 2 * spectral_height, plot_width, spectral_height]
        )
        # ax.yaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_major_locator(ticker.MultipleLocator(len(x_axis)))
        if np.any(bad_channels_orig):
            try:
                corrected_mean[bad_channels_orig] = 0
            except IndexError:
                pass
        ax.plot(x_axis, original_mean, lw=1, color="red")
        ax.plot(x_axis, corrected_mean, lw=1, color="green")
        ax.axhline(0, color="gray")
        ax.set_ylim(min_value, max_value)
        ax.text(
            0.03, 0.94, "corrected_mean", transform=ax.transAxes, ha="left", va="top"
        )
        for line_ in ax.get_xticklines() + ax.get_yticklines():
            line_.set_markersize(2)
        if i_drawn == 0:
            ax.set_title("PCA decomposition into Eigenspectra")
        i_drawn += 1

    for i, comp in enumerate(decomposition):
        ax = fig.add_axes(
            [
                padding,
                0.90 - ((i + 1 + extra_components) * spectral_height),
                plot_width,
                spectral_height,
            ]
        )
        if i_drawn == 0:
            ax.set_title("PCA decomposition into Eigenspectra")

        # ax.xaxis.set_major_locator(ticker.MultipleLocator(len(x_axis)))
        if i < number_components - 1:
            ax.xaxis.set_major_formatter(ticker.NullFormatter())
        else:
            ax.set_xlabel(x_axis_unit)
        color = "black"
        if i != 0:
            if decomposition_type == "pca":
                try:
                    if cutoff > pca.explained_variance_ratio_[i - 1]:
                        color = "gray"
                except AttributeError:
                    pass
        # if np.any(bad_channels):
        #     comp[bad_channels] = 0

        if channels_exclude:
            comp = np.concatenate(
                (
                    comp[0 : channels_exclude[0] - 1],
                    np.zeros(channels_exclude[1] - channels_exclude[0]),
                    comp[channels_exclude[0] - 1 :],
                )
            )
        ax.plot(x_axis, comp, "-", color=color, lw=1)

        if i == 0:
            label = "mean"
            explained_var = ""
            global_structure_var = ""
        else:
            label = "component %i" % i
            if decomposition_type == "pca":
                try:
                    explained_var = "{:1.3f}".format(
                        pca.explained_variance_ratio_[i - 1]
                    )
                except AttributeError:
                    pass
                global_structure_var = "{:1.3f}".format(global_structure[i - 1])
        ax.text(0.03, 0.94, label, transform=ax.transAxes, ha="left", va="top")
        ax.text(0.85, 0.94, explained_var, transform=ax.transAxes, ha="left", va="top")
        ax.text(
            0.5, 0.94, global_structure_var, transform=ax.transAxes, ha="left", va="top"
        )

        for line_ in ax.get_xticklines() + ax.get_yticklines():
            line_.set_markersize(2)
        i_drawn += 1
    return fig


def plot_example_correction(
    fig,
    example_correction,
    derived_scaled_noise_ratio,
    x_axis,
    config=None,
    plot_width=0.15,
    padding=0.05,
):
    x_axis = apply_pca_exclude_range(x_axis, config)
    fig.subplots_adjust(
        left=0.05, right=0.95, wspace=0.05, bottom=0.1, top=0.95, hspace=0.05
    )
    #for i in range(x_axis.shape[0]):
    #    print(x_axis[i], end=" ")
    print("x_axis shape: ", x_axis.shape)
    print (x_axis.shape[0])
    

    # Find the spectrum with the highest coefficient

    lowest_value = np.where(
        derived_scaled_noise_ratio == np.min(derived_scaled_noise_ratio)
    )
    try:
        lowest_value = lowest_value[0][0]
        first_example = example_correction[lowest_value]
    except KeyError as e:
        print(
            (
                "There is a problem with the example_correction dictionary. "
                "Key requested that does not exist."
            ),
            e,
        )
        raise SystemExit

    number_components = len(first_example["correction"])
    if number_components > 1:
        spectral_height = 0.8 / (number_components)
    else:
        spectral_height = 0.8

    for i in range(len(first_example["correction"])):
        if i not in first_example["correction"]:  # Check if key exists
            print(f"Warning: Missing correction for component {i}")
            continue
        correction = first_example["correction"][i]
        #print("size of y axis:")
        #print(len(correction["corrected_spec_scaled"]))
        t_ax = fig.add_axes(
            [
                5 * padding + 4 * plot_width,
                0.90 - (i + 1) * spectral_height,
                plot_width,
                spectral_height,
            ]
        )
        if i == 0:
            t_ax.set_title(
                "Example correction showing spectrum {}".format(int(lowest_value))
            )
        t_ax.plot(
            x_axis,
            copy.deepcopy(
                correction["original_spec_scaled"],
                config,
            ),
            "-k",
        )
        if correction["skipped"]:
            t_ax.text(
                0.03,
                0.94,
                (
                    "comp: {}; coeff: {:1.2f} not used, noise_ratio {:1.2f}; "
                    "noise_ratio_scaled {:1.2f}"
                ).format(
                    i + 1,
                    correction["coefficient"],
                    correction["noise_ratio"],
                    correction["noise_ratio_scaled"],
                ),
                transform=t_ax.transAxes,
                ha="left",
                va="top",
            )
            t_ax.plot(
                x_axis,
                copy.deepcopy(
                    correction["fitted_component_scaled"],
                    config,
                ),
                "-",
                color="red",
            )
            continue
        t_ax.plot(x_axis, correction["corrected_spec_scaled"], "-", color="gray")
        t_ax.plot(x_axis, correction["fitted_component_scaled"], color="blue")
        t_ax.text(
            0.03,
            0.94,
            (
                "comp: {}; coeff: {:1.2f}; "
                "noise_ratio {:1.2f}; noise_ratio_scaled {:1.2f}"
            ).format(
                i + 1,
                correction["coefficient"],
                correction["noise_ratio"],
                correction["noise_ratio_scaled"],
            ),
            transform=t_ax.transAxes,
            ha="left",
            va="top",
        )
        t_ax.set_xlabel("km/s")
    return fig


def plot_heatmap(
    fig,
    derived_scaled_noise_ratio,
    derived_scaled_noise_ratio_cut,
    derived_component_coefficients,
    cut_coefficients,
    x_axis,
    global_derived_scaled_noise_ratio,
    plot_width=0.15,
    padding=0.05,
):
    fig.subplots_adjust(
        left=0.05, right=0.95, wspace=0.05, bottom=0.1, top=0.95, hspace=0.05
    )
    ax = fig.add_axes([4 * padding + 3 * plot_width, 0.05, plot_width, 0.5])
    ax.set_title(
        "Coefficient values of components used in "
        "correction\nCutoff value {}".format(cut_coefficients)
    )
    # Minor ticks
    ax.set_xticks(
        np.arange(-0.5, len(derived_component_coefficients[0]), 1), minor=True
    )
    ax.set_yticks(
        np.arange(-0.5, len(derived_component_coefficients[0]), 1), minor=True
    )

    im1 = ax.imshow(
        1 / derived_scaled_noise_ratio,
        interpolation="none",
        aspect="auto",
        cmap="gray",
        alpha=0.5,
    )
    _ = ax.imshow(
        1 / derived_scaled_noise_ratio_cut,
        interpolation="none",
        aspect="auto",
        cmap="inferno",
        alpha=0.8,
    )
    # ax.hist(derived_scaled_noise_ratio[:, 0], bins=100)
    # , extent = [
    #                0, len(derived_component_coefficients[0]), 0, len(data)])
    ax.grid(which="minor", axis="x", color="w", linestyle="-", linewidth=2)
    ax.figure.colorbar(im1, ax=ax, pad=0)
    return fig


def gaussian(x, ampl, center, dev):
    """Computes the Gaussian function.

    Parameters
    ----------
    x : number
        Point to evaluate the Gaussian for.
    a : number
        Amplitude.
    b : number
        Center.
    c : number
        Width.

    Returns
    -------
    float
        Value of the specified Gaussian at *x*
    """
    eps = np.finfo(float).eps
    return ampl * np.exp(-((x - float(center)) ** 2) / (2.0 * dev**2 + eps))


def gaussian_fit(x, y, center_only=True, width=5):
    """Performs a Gaussian fitting of the specified data.

    Parameters
    ----------
    x : ndarray
        Data on the x axis.
    y : ndarray
        Data on the y axis.
    center_only: bool
        If True, returns only the center of the Gaussian
        for `interpolate` compatibility

    Returns
    -------
    ndarray or float
        If center_only is `False`, returns the parameters of the Gaussian
        that fits the specified data
        If center_only is `True`, returns the center position of the Gaussian
    """
    from scipy import optimize

    if len(x) < 3:
        # used RuntimeError to match errors raised in scipy.optimize
        raise RuntimeError("At least 3 points required for Gaussian fitting")

    initial = [np.max(y), x[0], (x[1] - x[0]) * width]
    params, pcov = optimize.curve_fit(gaussian, x, y, initial)
    if center_only:
        return params[1]
    else:
        return params


def plot_histogram(first_kde_gauss_peak_params, title=None):
    from matplotlib.figure import Figure

    fig = Figure(figsize=(25, 10))
    n_comp = len(first_kde_gauss_peak_params)

    for i in range(n_comp):
        ax1 = fig.add_subplot(
            math.ceil(math.sqrt(n_comp)), math.ceil(math.sqrt(n_comp)), i + 1
        )
        (
            amplitude,
            center,
            dev,
            X_plot,
            log_dens,
            indexes,
            gaussian_width,
        ) = first_kde_gauss_peak_params[i]

        ax1.plot(X_plot[:, 0], np.exp(log_dens), "-")

        if len(indexes) > 0:
            first_index = indexes[0]

            slice_ = slice(
                first_index - gaussian_width, first_index + gaussian_width + 1
            )
            try:
                ax1.plot(
                    X_plot[:, 0][slice_],
                    gaussian(X_plot[:, 0][slice_], amplitude, center, dev),
                    "-",
                )
                ax1.set_title(
                    "Scaled Noise ratio for component {}\n"
                    "{:1.2f} {:1.2f} {:1.2f}".format(i + 1, amplitude, center, dev)
                )
            except TypeError:
                print("ERror with gaus params", amplitude, center, dev)
                pass
        else:
            continue
        try:
            for ind in indexes:
                ax1.axvline(ind)
        except IndexError:
            pass
        # ax1.set_xscale('log')
    fig.savefig(title)


def plot_summary(
    original_spectra=None,
    corrected_spectra=None,
    pca_comp=None,
    sky_diff_array=None,
    pca=None,
    image_path=None,
    number_components=None,
    x_axis=None,
    x_axis_unit="km/s",
    cutoff=None,
    title=None,
    suptitle=None,
    channels_exclude=None,
    bad_value=None,
    good_channels_orig=None,
    bad_channels_orig=None,
    good_channels=None,
    bad_channels=None,
    norms=None,
    cut_coefficients=0.1,
    derived_component_coefficients=None,
    derived_scaled_noise_ratio=None,
    derived_scaled_noise_ratio_cut=None,
    example_correction=None,
    global_structure=None,
    global_derived_scaled_noise_ratio=None,
    plot_width=0.17,
    padding=0.02,
    contoured_plot=None,
    masked_plot=None,
    first_kde_gauss_peak_params=None,
    config=None,
):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=(25, 10))

    fig.suptitle(suptitle)

    canvas = FigureCanvasAgg(fig)  # This is intrinsically used by matplotlib

    fig = plot_pca_decomposition(
        fig=fig,
        decomposition=pca_comp,
        pca=pca,
        number_components=number_components,
        x_axis=x_axis,
        original_spectra=original_spectra,
        corrected_spectra=corrected_spectra,
        cutoff=cutoff,
        channels_exclude=channels_exclude,
        bad_channels_orig=bad_channels_orig,
        good_channels_orig=good_channels_orig,
        bad_channels=bad_channels,
        good_channels=good_channels,
        norms=norms,
        cut_coefficients=cut_coefficients,
        global_structure=global_structure,
        plot_width=plot_width,
        padding=padding,
    )

    ax1 = fig.add_axes([2 * padding + plot_width, 0.65, plot_width, 0.25])
    ax2 = fig.add_axes([2 * padding + plot_width, 0.05, plot_width, 0.5])
    ax3 = fig.add_axes([3 * padding + 2 * plot_width, 0.65, plot_width, 0.25])
    ax4 = fig.add_axes([3 * padding + 2 * plot_width, 0.05, plot_width, 0.5])
    ax5 = fig.add_axes([4 * padding + 3 * plot_width, 0.65, plot_width, 0.25])

    ax1_title = "Spectra used to derive components"

    im1 = ax1.imshow(
        sky_diff_array,
        aspect="auto",
        interpolation="none",
        extent=[x_axis[0], x_axis[-1], len(sky_diff_array), 0],
    )
    fig.colorbar(im1, ax=ax1, pad=0)
    ax1.set_title(ax1_title)
    ax1.set_xlabel(x_axis_unit)

    # Plot the corrected data
    ax2_title = "Spectra before PCA correction "
    if np.any(bad_channels):
        ax2_title += " WARNING: blanked channels present"
        vmin = ma.min(
            np.asarray(original_spectra)[
                np.where(~np.isclose(original_spectra, bad_value, atol=1e6))
            ]
        )
        vmax = ma.max(
            np.asarray(original_spectra)[
                np.where(~np.isclose(original_spectra, bad_value, atol=1e6))
            ]
        )
    else:
        vmin = ma.min(original_spectra)
        vmax = ma.max(original_spectra)
    # im2 = ax2.imshow(original_spectra, aspect="auto", interpolation='none',
    #                  extent=[x_axis[0], x_axis[-1], 0, len(original_spectra)])
    im2 = ax2.imshow(
        contoured_plot,
        aspect="auto",
        interpolation="none",
        extent=[x_axis[0], x_axis[-1], len(contoured_plot), 0],
    )
    fig.colorbar(im2, ax=ax2, pad=0)
    ax2.set_title(ax2_title)
    ax2.set_xlabel(x_axis_unit)

    ax3_title = "Difference Original vs Corrected"
    if np.any(bad_channels):
        ax3_title += " WARNING: blanked channels present"
        vmin = ma.min(
            corrected_spectra[
                np.where(~np.isclose(corrected_spectra, bad_value, atol=1e6))
            ]
        )
        vmax = ma.max(
            corrected_spectra[
                np.where(~np.isclose(corrected_spectra, bad_value, atol=1e6))
            ]
        )
    else:
        vmin = ma.min(corrected_spectra)
        vmax = ma.max(corrected_spectra)

    im3 = ax3.imshow(
        original_spectra - corrected_spectra,
        aspect="auto",
        interpolation="none",
        extent=[
            x_axis[0],
            x_axis[-1],
            len(corrected_spectra),
            0,
        ],
    )
    fig.colorbar(im3, ax=ax3, pad=0)
    ax3.set_title(ax3_title)
    ax3.set_xlabel(x_axis_unit)

    ax4_title = "Spectra after PCA correction"
    im4 = ax4.imshow(
        corrected_spectra,
        aspect="auto",
        interpolation="none",
        extent=[
            x_axis[0],
            x_axis[-1],
            len(corrected_spectra),
            0,
        ],
    )

    ax4.set_title(ax4_title)
    ax4.set_xlabel(x_axis_unit)
    fig.colorbar(im4, ax=ax4, pad=0)

    ax5_title = "Masked array used to fit components to spectrum"
    _ = ax5.imshow(
        masked_plot,
        aspect="auto",
        interpolation="none",
        extent=[x_axis[0], x_axis[-1], len(masked_plot), 0],
    )
    ax5.set_title(ax5_title)
    ax5.set_xlabel(x_axis_unit)

    fig = plot_example_correction(
        fig,
        example_correction,
        derived_scaled_noise_ratio,
        x_axis,
        config,
        plot_width=plot_width,
        padding=padding,
    )

    fig = plot_heatmap(
        fig,
        derived_component_coefficients=derived_component_coefficients,
        derived_scaled_noise_ratio=derived_scaled_noise_ratio,
        derived_scaled_noise_ratio_cut=derived_scaled_noise_ratio_cut,
        cut_coefficients=cut_coefficients,
        x_axis=x_axis,
        global_derived_scaled_noise_ratio=global_derived_scaled_noise_ratio,
        plot_width=plot_width,
        padding=padding,
    )

    plot_histogram(first_kde_gauss_peak_params, image_path.replace(".png", "_2.png"))

    fig.savefig(image_path)


def command_line_arguments():
    parser = OptionParser()

    parser = add_correct_options(parser)
    parser = add_common_options(parser)

    try:
        (options, args) = parser.parse_args()
    except KeyError:
        pyclass.message(pyclass.seve.e, "PCA", "Invalid option")
        pyclass.sicerror()
        raise InvalidOption

    # Start with command line options
    config = options.__dict__
    print(config)

    # Setup parameters
    if options.config_file:
        pyclass.message(
            pyclass.seve.i,
            "PCA",
            "Read configuration from {}".format(options.config_file),
        )
        config_from_file = load_config_file(options.config_file)
        if not config_from_file:
            raise InvalidOption
        if options.debug:
            print(config_from_file)

        # Update with config file settings, but only if not set in command line
        if "common" in config_from_file:
            for key, value in config_from_file["common"].items():
        #        if key not in config:
                config[key] = value
        #        else:
        #            pyclass.message(
        #                pyclass.seve.i,
        #                "PCA",
        #                "Command line option takes precedence over config file for: {}".format(
        #                    key
        #                ),
        #            )
#
        if "correct" in config_from_file:
            for key, value in config_from_file["correct"].items(): 
        #       if key not in config:
                config[key] = value
        #        else:
        #            pyclass.message(
        #                pyclass.seve.i,
        #                "PCA",
        #                "Command line option takes precedence over config file for: {}".format(
        #                    key
        #                ),
        #            )

    if options.input_file is None:
        pyclass.message(
            pyclass.seve.e, "PCA", "The input file name is needed: /input_file FILENAME"
        )
        pyclass.sicerror()
        raise MissingMandatoryOption

    if not pyclass.gotgdict():
        pyclass.get(verbose=False)

    try:
        pyclass.comm("def real noise_measure")
    except pgutils.PygildasError:
        pass

    config["line_windows"] = get_line_windows(config)
    if config["output_folder"] and not os.path.isdir(config["output_folder"]):
        subprocess.call(["mkdir", "-p", config["output_folder"]])

    return config


def prepare_spectrum_for_plot(spectrum, config, scaler=None):
    channels_exclude = get_exclude_channels(config)
    return_spectrum = copy.deepcopy(spectrum)
    if scaler is not None:
        return_spectrum = scaler.inverse_transform(return_spectrum.reshape(1, -1))[0]
    if channels_exclude:
        return_spectrum = np.concatenate(
            (
                return_spectrum[0 : channels_exclude[0] - 1],
                np.zeros(channels_exclude[1] - channels_exclude[0]),
                return_spectrum[channels_exclude[0] - 1 :],
            )
        )
    #for i in range(len(return_spectrum)):
    #        print(return_spectrum[i])
    return return_spectrum


def derive_noise_ratio_cutoff(derived_scaled_noise_ratio):
    import peakutils
    from sklearn.impute import SimpleImputer
    from sklearn.neighbors import KernelDensity

    first_kde_gauss_peak_params = {}

    for i, comp in enumerate(derived_scaled_noise_ratio.T):
        comp = np.abs(comp)
        X = comp[:, np.newaxis]
        X_plot = np.linspace(0, 60, 60)[:, np.newaxis]
        kde = KernelDensity(kernel="gaussian", bandwidth=1.5).fit(X)
        log_dens = kde.score_samples(X_plot)

        indexes = peakutils.indexes(np.exp(log_dens))
        if len(indexes) > 0:
            gaussian_width = 2
            first_index = indexes[0]
            first_slice = first_index - gaussian_width
            if (first_slice) < 0:
                first_slice = 0
            slice_ = slice(first_slice, first_index + gaussian_width + 1)
            try:
                amplitude, center, dev = gaussian_fit(
                    X_plot[:, 0][slice_],
                    np.exp(log_dens[slice_]),
                    center_only=False,
                    width=gaussian_width,
                )
            except RuntimeError as e:
                print(e)
                amplitude, center, dev = (None, None, None)
        else:
            gaussian_width = None
            amplitude, center, dev = (None, None, None)
        first_kde_gauss_peak_params[i] = [
            amplitude,
            center,
            dev,
            X_plot,
            log_dens,
            indexes,
            gaussian_width,
        ]
    return first_kde_gauss_peak_params


def should_component_be_used(gauss_params, config):
    """Calculates a flag if this component should be used in general

    for correction and the noise_ratio cutoff value above which a
    component is not used.
    """
    amplitude, center, dev, X_plot, log_dens, indexes, width = gauss_params
    use_component = True
    global_cutoff = False
    if center and (int(config["global_noise_ratio_cutoff"]) != 0):
        use_component = center < float(config["global_noise_ratio_cutoff"])
        global_cutoff = True
    else:
        use_component = True
    if use_component and global_cutoff:
        # If there is a global cutoff defined and the component should be used
        # accecpt all spectra that are below the cutoff of 3 sigma away from
        # the gaussion of the first peak of the distribution of noise_ratios
        cutoff_value = center + 3 * dev
    else:
        cutoff_value = None
    return use_component, cutoff_value


def derive_noise_ratio_and_windows(
    config,
    threshold=None,
    numbers=None,
    scaler=None,
    spectrum_collection=None,
    decomposition=None,
    subscan=None,
):
    """This function"""

    prelim_corrected_spectra = None
    derived_scaled_noise_ratio = np.zeros(
        [len(spectrum_collection["numbers"]), len(decomposition.pca.components_)]
    )

    this_components = np.asarray(decomposition.pca.components_)
    for i, number in enumerate(numbers):
        if threshold is None:
            this_threshold = threshold
        elif len(threshold[i, :]) == 1:
            this_threshold = threshold
        else:
            this_threshold = threshold[i, :]
        if subscan:
            this_subscan = spectrum_collection[number]["subscan"]
            if this_subscan != subscan:
                continue
        bad_channels_mask = spectrum_collection[number]["bad_channels_mask"]
        bad_channels = spectrum_collection[number]["bad_channels"]
        good_channels = spectrum_collection[number]["good_channels"]
        spectrum = copy.deepcopy(spectrum_collection[number]["exclude_range"])

        if config["do_scale"] and scaler:
            new_spec = scaler.transform(spectrum.reshape(1, -1))[0]
            scaled_fit_spec = scaler.transform(spectrum.reshape(1, -1))[0]
        else:
            new_spec = copy.deepcopy(spectrum)
            scaled_fit_spec = copy.deepcopy(spectrum)

        if this_threshold is not None:
            scaled_fit_spec = ma.array(scaled_fit_spec, mask=this_threshold, copy=True)
        spectrum_collection[number]["new_spec"] = copy.deepcopy(new_spec)
        spectrum_collection[number]["scaled_fit_spec"] = copy.deepcopy(scaled_fit_spec)

        dot_calculation_failed = False
        if np.any(bad_channels):
            try:
                coeff = ma.dot(
                    this_components[:, spectrum_collection[number]["good_channels"]],
                    scaled_fit_spec[spectrum_collection[number]["good_channels"]],
                )
            except ValueError as e:
                print(e)
                dot_calculation_failed = True
        else:
            try:
                coeff = ma.dot(this_components, scaled_fit_spec)
            except ValueError as e:
                print(e)
                dot_calculation_failed = True
        if dot_calculation_failed:
            print(
                "It looks like you are trying to use components that have "
                "a different shape than your spectra. Did you change the "
                "requested dimensions between creating the components and "
                "executing the corrections? You can try to remove the hidden "
                ".pca* files to start from scratch and rerun decomposition and "
                "correction."
            )
            raise SystemExit
        spectrum_collection[number]["coeff"] = copy.deepcopy(coeff)
        spectrum_collection[number]["components"] = {}
        for j, component in enumerate(this_components):
            spectrum_collection[number]["components"][j] = {}

            # Scale this component to the level of the spectrum
            scaled_component = copy.deepcopy((np.dot(coeff[j], component)))

            # Store the scaled component
            spectrum_collection[number]["components"][j]["scaled_component"] = (
                copy.deepcopy(scaled_component)
            )

            noise_in_scaled_component = scaled_component.std()

            if this_threshold is not None:
                masked_new_spec = ma.masked_array(new_spec, this_threshold, copy=True)
            else:
                masked_new_spec = new_spec

            # To smooth out high resolution fluctuations and create a smooth
            # component not the data is smoothed if "smoothing_kernel_size"
            # is configured. Here the spectrum is smoothed to be able to
            # compare this to the noise in the smoothed componentn
            if config["smoothing_kernel_size"]:
                kernel_size = config["smoothing_kernel_size"]
                kernel = np.ones(kernel_size) / kernel_size
                new_smoothed_spec = np.convolve(new_spec, kernel, mode="same")
            else:
                new_smoothed_spec = new_spec

            if this_threshold is not None:
                masked_new_smoothed_spec = ma.masked_array(
                    new_smoothed_spec, this_threshold, copy=True
                )
            else:
                masked_new_smoothed_spec = new_smoothed_spec

            noise_in_spectrum_smoothed = ma.std(masked_new_smoothed_spec)

            noise_in_spectrum = ma.std(masked_new_spec)

            noise_ratio = noise_in_spectrum / noise_in_scaled_component

            noise_ratio_scaled = noise_in_spectrum_smoothed / noise_in_scaled_component

            derived_scaled_noise_ratio[i][j] = copy.deepcopy(noise_ratio_scaled)

            spectrum_collection[number]["components"][j]["noise_ratio_scaled"] = (
                copy.deepcopy(noise_ratio_scaled)
            )
            spectrum_collection[number]["components"][j]["noise_ratio"] = noise_ratio
            new_spec = np.subtract(
                new_spec, scaled_component, where=good_channels, out=new_spec
            )
        new_spec = ma.masked_array(new_spec, bad_channels_mask, copy=True)
        if prelim_corrected_spectra is None:
            prelim_corrected_spectra = new_spec
        else:
            prelim_corrected_spectra = ma.vstack((prelim_corrected_spectra, new_spec))
    return (prelim_corrected_spectra, spectrum_collection, derived_scaled_noise_ratio)


def pca_correct(config=None):
    if config is None:
        try:
            config = command_line_arguments()
        except InvalidOption:
            return
        except MissingMandatoryOption:
            return
    # Input parameters

    df = create_index(config)

    set_line_windows(config)
    component_map = []

    # The pca reduction has to tbe done per frontend, because each of the
    # frontends has a different response
    groups = df.groupby(["scan", "telescope", "mission_id"])

    # Iterate over individual groups
    for name, scan_group in groups:
        scan = str(name[0]).strip()

        # Only work on scans that are requested or all
        if config["limit_to_scan"]:
            if int(config["limit_to_scan"]) != int(scan):
                continue

        telescope = str(name[1]).strip()
        # Only work on telescopes that are requested or all
        if config["limit_to_telescope"]:
            if config["limit_to_telescope"] != telescope:
                continue

        mission_id = str(name[2]).strip()

        config["telescope"] = telescope
        config["mission_id"] = mission_id
        if config["no_subscan_grouping"]:
            subscan_groups = scan_group
        else:
            subscan_groups = scan_group.groupby(["subscan"])

        original_specs = []
        corrected_specs = None
        final_contoured_plot = None
        final_derived_scaled_noise_ratio = None
        final_derived_scaled_noise_ratio_cut = None
        final_masked_spec_array_to_correct = None
        example_correction = {}
        scan_iteration = -1

        for name, group in subscan_groups:
            if config["no_subscan_grouping"]:
                subscan = "1"
            else:
                subscan = str(name).strip()
            # Find the name for the component file for this group
            # these files are created per flight and telescope
            try:
                science_source = get_and_check_science_sources(
                    group, scan, subscan, config
                )
            except NoScienceSourceFound:
                continue
            except MultipleScienceSourcesFound:
                continue

            decomposition_pickle_file = create_pickle_file_name(config)

            # Timer to measure time taken
            start = time.time()
            try:
                with open(decomposition_pickle_file, "rb") as filehandler:
                    decomposition = pickle.load(filehandler)
            except IOError:
                print("Loading pickled index failed: ", decomposition_pickle_file)
                continue

            print("Loading pickled index took: ", time.time() - start)
            pyclass.message(
                pyclass.seve.i,
                "PCA",
                "Processing scan: {}; subscan: {}; pixel: {}".format(
                    scan, subscan, telescope
                ),
            )

            numbers = group.loc[group["source"] == science_source, "number"].tolist()

            pyclass.message(
                pyclass.seve.i,
                "PCA",
                "Removing components that explain more than "
                "{:1.0f}% of the variance in the {} data".format(
                    config["cutoff"] * 100, config["pca_source"]
                ),
            )
            pyclass.message(pyclass.seve.i, "PCA", "Version 1.1")

            # Dereive the mean spectrum for this group
            spec_array_to_correct = None
            bad_channel_array = None
            good_channel_array = None
            spectrum_collection = {}
            spectrum_collection["numbers"] = numbers
            bad_channel_mask_array = None

            # 1. Load the all spectra for this group into an array to
            # use for line_detection
            for number in numbers:
                # subscan = group.loc[group["number"] == number].subscan.values
                # subscan = int(subscan)
                # if subscan not in subscan_arrays:
                #     subscan_arrays[subscan] = None
                # if subscan not in subscans:
                #     subscans += [subscan]
                spectrum_collection[number] = {}
                try:
                    spectrum, outlier_mask = prepare_spectrum(number, config)
                except pgutils.PygildasError:
                    continue
                spectrum_rms = copy.deepcopy(pyclass.gdict.sigma.__sicdata__)

                spectrum_collection[number]["original_rms"] = copy.deepcopy(
                    spectrum_rms
                )
                spectrum_collection[number]["original"] = copy.deepcopy(spectrum)

                bad_value = float(
                    copy.deepcopy(pyclass.gdict.r.head.spe.bad.__sicdata__)
                )
                spectrum_collection[number]["despiked"] = copy.deepcopy(spectrum)

                bad_channels_mask = np.where(
                    np.isclose(spectrum, bad_value, atol=1e6), 1, 0
                )
                if bad_channel_mask_array is None:
                    bad_channel_mask_array = [bad_channels_mask]
                else:
                    bad_channel_mask_array += [bad_channels_mask]

                original_specs += [
                    ma.masked_array(spectrum, bad_channels_mask, copy=True)
                ]

                bad_channel_mask_array = None
                spectrum = apply_pca_exclude_range(spectrum, config)

                bad_channels_mask = np.where(
                    np.isclose(spectrum, bad_value, atol=1e6), 1, 0
                )
                spectrum_collection[number]["bad_channels_mask"] = copy.deepcopy(
                    bad_channels_mask
                )
                if bad_channel_mask_array is None:
                    bad_channel_mask_array = [bad_channels_mask]
                else:
                    bad_channel_mask_array += [bad_channels_mask]

                # Now remove the excluded range from the spectrum
                # This makes the spectrum shorter.
                spectrum_collection[number]["exclude_range"] = copy.deepcopy(spectrum)

                if spec_array_to_correct is None:
                    spec_array_to_correct = ma.masked_array(
                        spectrum, bad_channels_mask, copy=True
                    )
                else:
                    spec_array_to_correct = ma.vstack(
                        (
                            spec_array_to_correct,
                            ma.masked_array(spectrum, bad_channels_mask, copy=True),
                        )
                    )

                bad_channels_orig = np.where(
                    np.isclose(spectrum, bad_value, atol=1e6), True, False
                )
                bad_channels = copy.deepcopy(bad_channels_orig)
                if bad_channel_array is None:
                    bad_channel_array = [bad_channels_orig]
                else:
                    bad_channel_array += [bad_channels_orig]

                good_channels_orig = ~bad_channels
                good_channels = copy.deepcopy(good_channels_orig)
                if good_channel_array is None:
                    good_channel_array = [bad_channels_orig]
                else:
                    good_channel_array += [bad_channels_orig]

                spectrum_collection[number]["good_channels"] = copy.deepcopy(
                    good_channels
                )
                spectrum_collection[number]["bad_channels"] = copy.deepcopy(
                    bad_channels
                )
                spectrum_collection[number]["reference"] = copy.deepcopy(
                    pyclass.gdict.reference.__sicdata__
                )
                spectrum_collection[number]["velocity"] = copy.deepcopy(
                    pyclass.gdict.velocity.__sicdata__
                )
                spectrum_collection[number]["velo_res"] = copy.deepcopy(
                    pyclass.gdict.velo_step.__sicdata__
                )
                spectrum_collection[number]["subscan"] = subscan

            # and the scaler for the working data
            print(config["do_scale"], type(config["do_scale"]))
            if config["do_scale"]:
                if len(spec_array_to_correct.shape) > 1:
                    scaler = preprocessing.StandardScaler().fit(spec_array_to_correct)
                else:
                    scaler = None
            else:
                scaler = None
            try:
                this_components = decomposition.pca.components_
            except AttributeError:
                print("NO Components")
                continue

            # Detect line windows

            spec_array_to_correct_original = copy.deepcopy(spec_array_to_correct)

            line_kernel_size = 51
            if config["line_kernel_size"]:
                line_kernel_size = int(config["line_kernel_size"])

            # threshold, contoured_plot = find_lines(
            #     spec_array_to_correct,
            #     config=config,
            #     kernel_size=line_kernel_size
            # )

            try:
                (
                    prelim_corrected_spectra,
                    spectrum_collection,
                    derived_scaled_noise_ratio,
                ) = derive_noise_ratio_and_windows(
                    config,
                    numbers=numbers,
                    scaler=scaler,
                    spectrum_collection=spectrum_collection,
                    decomposition=decomposition,
                    subscan=subscan,
                )
            except RuntimeWarning as e:
                print(e)
                print(
                    "This looks as if there is an empty slice "
                    "in the data. Trying to continue"
                )
                continue
            # if a subscan has only a single spectrum the line detection
            # algorithm and also the derived windows are in the data wrong format
            if len(np.shape(prelim_corrected_spectra)) == 1:
                prelim_corrected_spectra_ = []
                for _ in range(10):
                    prelim_corrected_spectra_.append(prelim_corrected_spectra)
                prelim_corrected_spectra = prelim_corrected_spectra_
            threshold, contoured_plot = find_lines(
                prelim_corrected_spectra,
                plot_image=spec_array_to_correct,
                kernel_size=line_kernel_size,
                config=config,
            )
            this_cutoff_std = 2
            this_line_kernel_size = line_kernel_size
            for iteration in [0, 1, 2]:
                if iteration == 2:
                    this_cutoff_std = 3
                # this_line_kernel_size = this_line_kernel_size - 10

                spec_array_to_correct_original_temp = copy.deepcopy(
                    spec_array_to_correct_original
                )
                (
                    prelim_corrected_spectra,
                    spectrum_collection,
                    derived_scaled_noise_ratio,
                ) = derive_noise_ratio_and_windows(
                    config,
                    numbers=numbers,
                    scaler=scaler,
                    threshold=threshold,
                    spectrum_collection=spectrum_collection,
                    decomposition=decomposition,
                    subscan=subscan,
                )
                if len(np.shape(prelim_corrected_spectra)) == 1:
                    prelim_corrected_spectra_ = []
                    for _ in range(10):
                        prelim_corrected_spectra_.append(prelim_corrected_spectra)
                    prelim_corrected_spectra = prelim_corrected_spectra_

                threshold, contoured_plot = find_lines(
                    prelim_corrected_spectra,
                    plot_image=spec_array_to_correct_original_temp,
                    kernel_size=this_line_kernel_size,
                    config=config,
                    cutoff_std=this_cutoff_std,
                )

            if final_contoured_plot is None:
                final_contoured_plot = contoured_plot
            else:
                final_contoured_plot = ma.vstack((final_contoured_plot, contoured_plot))
            if config["line_windows"]:
                flat_line_windows = np.asarray(config["line_windows"]).flatten()

                window_channels = [
                    derive_channel_from_velocity(
                        reference=spectrum_collection[number]["reference"],
                        central_velocity=spectrum_collection[number]["velocity"],
                        velo_res=spectrum_collection[number]["velo_res"],
                        channel_velocity=window_entry,
                    )
                    for window_entry in flat_line_windows
                ]

                window_channels = sorted(window_channels)
                window_channels = [0] + window_channels + [-1]
                for k in range(0, len(window_channels), 2):
                    threshold[:, window_channels[k] : window_channels[k + 1]] = 0

            # contoured_plot = ma.masked_array(
            #     contoured_plot, bad_channel_mask_array)

            global_derived_scaled_noise_ratio = np.median(
                derived_scaled_noise_ratio, axis=0
            )
            derived_component_coefficients = np.zeros(
                [len(numbers), len(decomposition.pca.components_)]
            )

            # Derive Kernel density estimation of the coefficient distribution
            # and fit a gaussian to the first peak
            first_kde_gauss_peak_params = derive_noise_ratio_cutoff(
                derived_scaled_noise_ratio
            )

            # Copy the scaled_noise_ratio array to allow to mark positions that
            # have been dropped based on the filtering
            derived_scaled_noise_ratio_cut = copy.deepcopy(derived_scaled_noise_ratio)

            # Dictionary to collect the reduction steps to show what the PCA
            # does
            masked_spec_array_to_correct = None
            for i, number in enumerate(numbers):
                scan_iteration += 1
                if len(threshold[i]) <= 1:
                    this_threshold = copy.deepcopy(threshold)
                else:
                    this_threshold = copy.deepcopy(threshold[i])
                # Derive the channels that correspond to the line
                # velocities
                example_correction[scan_iteration] = {}
                example_correction[scan_iteration]["length"] = len(
                    decomposition.pca.components_
                )
                example_correction[scan_iteration]["correction"] = {}

                spectrum = spectrum_collection[number]["original"]

                # Save an orignal copy of the spectrum as a Reference
                # to e.g. later fill the excluded channels back in
                new_spec_orig = copy.deepcopy(spectrum)

                # Now apply the exclude range if defined. From here on
                # our spectra are shorter
                bad_channels_mask = spectrum_collection[number]["bad_channels_mask"]
                bad_channels = spectrum_collection[number]["bad_channels"]
                good_channels = spectrum_collection[number]["good_channels"]
                spectrum = spectrum_collection[number]["exclude_range"]
                if config["do_scale"]:
                    if scaler is not None:
                        new_spec = scaler.transform(spectrum.reshape(1, -1))[0]
                        scaled_fit_spec = scaler.transform(spectrum.reshape(1, -1))[0]
                    else:
                        new_spec = spectrum
                        scaled_fit_spec = spectrum
                else:
                    new_spec = spectrum
                    scaled_fit_spec = spectrum
                apply_mask = np.asarray(bad_channels_mask, dtype=np.uint8)
                if threshold is not None:
                    apply_mask = ma.mask_or(this_threshold, apply_mask)

                # apply_mask = threshold[i]
                scaled_fit_spec = ma.array(scaled_fit_spec, mask=apply_mask, copy=True)

                # To account for strong lines that might disturb the
                # adjustment of the pca components to the spectrum
                # we fill emission lines with artificial noise if this
                # is requested on the command line
                if masked_spec_array_to_correct is None:
                    masked_spec_array_to_correct = scaled_fit_spec
                else:
                    masked_spec_array_to_correct = ma.vstack(
                        (masked_spec_array_to_correct, scaled_fit_spec)
                    )

                # Read the value that is considered to mark blanked channels
                channels_exclude = get_exclude_channels(config)

                coeff = spectrum_collection[number]["coeff"]
                # Save the coeffiecients for later usage (e.g. plotting)
                derived_component_coefficients[i] = coeff

                # Load the components for these the channels have
                # already been excluded
                components_used_for_correction = []
                components_used_for_correction_list = []

                for j, component in enumerate(this_components):
                    component_dict = spectrum_collection[number]["components"][j]
                    scaled_component = component_dict["scaled_component"]
                    noise_ratio = component_dict["noise_ratio"]
                    noise_ratio_scaled = component_dict["noise_ratio_scaled"]

                    use_component, cutoff_value = should_component_be_used(
                        first_kde_gauss_peak_params[j], config
                    )

                    # Noise ratio cutoff is Boolean flag that determines if
                    # the component is ignored based o nthe noise_ratio
                    noise_ratio_cutoff = True
                    if use_component:
                        if int(config["global_noise_ratio_cutoff"]) != 0:
                            if config["noise_ratio_cutoff"] == 0 and cutoff_value:
                                noise_ratio_cutoff = noise_ratio_scaled > cutoff_value
                            else:
                                #
                                noise_ratio_cutoff = (
                                    noise_ratio_scaled > config["noise_ratio_cutoff"]
                                )
                        else:
                            noise_ratio_cutoff = False
                    coefficient_cutoff = False
                    if float(config["cut_coefficients"]) != 0:
                        coefficient_cutoff = abs(coeff[j]) < float(
                            config["cut_coefficients"]
                        )

                    if coefficient_cutoff or noise_ratio_cutoff:
                        derived_scaled_noise_ratio_cut[i][j] = np.nan
                        print("Skip component {} to small-> continue".format(j))
                        plot_scaled_component = prepare_spectrum_for_plot(
                            scaled_component, config, scaler
                        )
                        plot_corrected_spec = prepare_spectrum_for_plot(
                            new_spec, config, scaler
                        )
                        plot_original_spec = prepare_spectrum_for_plot(
                            new_spec, config, scaler
                        )
                        example_correction[scan_iteration]["correction"][j] = {
                            "skipped": True,
                            "coefficient": coeff[j],
                            "fitted_component": plot_scaled_component,
                            "corrected_spec": plot_corrected_spec,
                            "original_spec": plot_original_spec,
                            "corrected_spec_scaled": new_spec,
                            "original_spec_scaled": new_spec,
                            "fitted_component_scaled": scaled_component,
                            "noise_ratio": noise_ratio,
                            "noise_ratio_scaled": noise_ratio_scaled,
                            "full_original": spectrum_collection[number]["original"],
                        }
                    else:
                        components_used_for_correction += [[j, component]]
                        components_used_for_correction_list += [component]
                        print("component used: {}".format(abs(coeff[j])))

                # If there are no components found that should be applied
                # we write back the original spectrum

                if not components_used_for_correction:
                    if config["do_scale"] and scaler:
                        this_new_spec = scaler.inverse_transform(
                            new_spec.reshape(1, -1)
                        )[0]
                    else:
                        this_new_spec = new_spec
                    # new_spec = new_spec * spectrum_rms
                    # and fill back in the excluded channels unchanged
                    this_new_spec_final = refill_pca_exclude_range(
                        this_new_spec, new_spec_orig, config
                    )
                    # Write the changed spectrum back to the ry variable to
                    # save the corrected version of the spectrum
                    spectrum, outlier_mask = prepare_spectrum(
                        number,
                        config,
                        replace_ry=this_new_spec_final,
                        z_score_blank=10,
                        baseline_order=0,
                    )
                    pyclass.comm("write")
                    bad_value = float(
                        copy.deepcopy(pyclass.gdict.r.head.spe.bad.__sicdata__)
                    )

                    bad_channels_mask = np.where(
                        np.isclose(spectrum, bad_value, atol=1e6), 1, 0
                    )
                    spectrum = ma.masked_array(spectrum, bad_channels_mask, copy=True)

                    if corrected_specs is None:
                        corrected_specs = spectrum
                    else:
                        corrected_specs = ma.vstack((corrected_specs, spectrum))
                    example_correction[scan_iteration]["final"] = {
                        "final_spec": copy.deepcopy(spectrum),
                        "scaled_final_spec": copy.deepcopy(new_spec),
                    }
                    continue

                # If we got here there are components to be used for
                # corrections to make the best fit to the data we recalculate
                # the coefficients for the subset of the components that are
                # used for this spectrum
                if np.any(bad_channels):
                    components_used_for_correction_list = np.asarray(
                        components_used_for_correction_list
                    )
                    coeff = ma.dot(
                        components_used_for_correction_list[:, good_channels],
                        scaled_fit_spec[good_channels],
                    )
                else:
                    coeff = ma.dot(components_used_for_correction_list, scaled_fit_spec)

                # Now we aplly the adjuxsted component
                for j, component in enumerate(components_used_for_correction):
                    this_component_number, component = component
                    component_dict = spectrum_collection[number]["components"][
                        this_component_number
                    ]
                    last_new_spec = copy.deepcopy(new_spec)

                    scaled_component = copy.deepcopy((ma.dot(coeff[j], component)))

                    noise_in_component = component.std()
                    noise_in_scaled_component = scaled_component.std()
                    noise_in_spectrum = new_spec.std()
                    noise_ratio = noise_in_spectrum / noise_in_component
                    noise_ratio_scaled = noise_in_spectrum / noise_in_scaled_component

                    noise_ratio = component_dict["noise_ratio"]
                    noise_ratio_scaled = component_dict["noise_ratio_scaled"]

                    if config["cutoff"]:
                        if (
                            decomposition.pca.explained_variance_ratio_[j]
                            <= config["cutoff"]
                        ):
                            continue

                    # Only subtract channels that are not blanked
                    new_spec = np.subtract(
                        new_spec, scaled_component, where=good_channels, out=new_spec
                    )
                    # new_spec_mean = new_spec[good_channels].mean(0)
                    # new_spec = np.subtract(new_spec, new_spec_mean,
                    #                        where=good_channels,
                    #                        out=new_spec)

                    plot_scaled_component = prepare_spectrum_for_plot(
                        scaled_component, config, scaler
                    )
                    plot_corrected_spec = prepare_spectrum_for_plot(
                        new_spec,
                        config,
                        scaler,
                    )
                    plot_original_spec = prepare_spectrum_for_plot(
                        last_new_spec,
                        config,
                        scaler,
                    )
                    example_correction[scan_iteration]["correction"][
                        this_component_number
                    ] = {
                        "skipped": False,
                        "coefficient": coeff[j],
                        "fitted_component": plot_scaled_component,
                        "corrected_spec": plot_corrected_spec,
                        "original_spec": plot_original_spec,
                        "corrected_spec_scaled": new_spec,
                        "original_spec_scaled": last_new_spec,
                        "fitted_component_scaled": scaled_component,
                        "noise_ratio": noise_ratio,
                        "noise_ratio_scaled": noise_ratio_scaled,
                    }

                # if config["custom_baseline_fit"]:
                #     pyclass.message(
                #         pyclass.seve.i,
                #         "PCA",
                #         "Perfoming cusomg baseline fit of order {}".format(
                #             config["custom_baseline_fit"])
                #     )
                #     final_fit_spec = copy.deepcopy(new_spec)

                #     final_fit_spec = ma.array(final_fit_spec,
                #                               mask=this_threshold)

                #     xp = np.linspace(
                #         0,
                #         len(final_fit_spec), len(final_fit_spec))
                #     z = np.poly1d(
                #         ma.polyfit(
                #             xp,
                #             final_fit_spec,
                #             int(config["custom_baseline_fit"]))
                #         )
                #     new_spec = np.subtract(new_spec, z(xp),
                #                                     where=good_channels,
                #                                     out=new_spec)

                # mean = ma.mean(new_spec)
                # new_spec -= mean

                # example_correction[i]["correction"][
                #     len(example_correction[i]["correction"])] = {
                #     "skipped": False,
                #     "coefficient": coeff[j],
                #     "fitted_component": z(xp),
                #     "corrected_spec": plot_corrected_spec,
                #     "original_spec": plot_original_spec,
                #     "corrected_spec_scaled": new_spec,
                #     "original_spec_scaled": last_new_spec,
                #     "fitted_component_scaled": z(xp),
                #     "noise_ratio": noise_ratio,
                #     "noise_ratio_scaled": noise_ratio_scaled
                # }

                # Now scale the data back to its original scale
                if config["do_scale"] and scaler:
                    new_spec = scaler.inverse_transform(new_spec.reshape(1, -1))[0]
                # new_spec = new_spec * spectrum_rms
                # and fill back in the excluded channels unchanged
                new_spec_final = refill_pca_exclude_range(
                    new_spec, new_spec_orig, config
                )
                # Write the changed spectrum back to the ry variable to
                # save the corrected version of the spectrum
                baseline_order = 0
                if config["custom_baseline_fit"]:
                    baseline_order = int(config["custom_baseline_fit"])

                threshold_default = copy.deepcopy(new_spec_orig)
                threshold_default[:] = 0
                this_threshold = refill_pca_exclude_range(
                    this_threshold, threshold_default, config
                )

                spectrum, outlier_mask = prepare_spectrum(
                    number,
                    config=config,
                    replace_ry=new_spec_final,
                    window_mask=this_threshold,
                    baseline_order=baseline_order,
                    z_score_blank=10,
                )

                pyclass.comm("write")

                if config["component_map"]:
                    this_component_map = {}
                    this_component_map["lambda"] = copy.deepcopy(
                        pyclass.gdict.r.head.pos.lam.__sicdata__
                    )
                    this_component_map["beta"] = copy.deepcopy(
                        pyclass.gdict.r.head.pos.lam.__sicdata__
                    )

                    for x, _ in enumerate(this_components):
                        this_component_map[x] = coeff[x]
                    component_map.append(this_component_map)
                if config["create_self_pca"]:
                    if np.sum(this_threshold) == 0:
                        pyclass.comm("mod source SELF-PCA")
                        pyclass.comm("write")

                bad_value = float(
                    copy.deepcopy(pyclass.gdict.r.head.spe.bad.__sicdata__)
                )

                bad_channels_mask = np.where(
                    np.isclose(spectrum, bad_value, atol=1e6), 1, 0
                )
                spectrum = ma.masked_array(spectrum, bad_channels_mask, copy=True)
                if corrected_specs is None:
                    corrected_specs = spectrum
                else:
                    corrected_specs = ma.vstack((corrected_specs, spectrum))

                example_correction[scan_iteration]["final"] = {
                    "final_spec": copy.deepcopy(spectrum),
                    "scaled_final_spec": copy.deepcopy(new_spec),
                }
                # scan_iteration += 1
                # Also write out the components to the class file

            if final_derived_scaled_noise_ratio is None:
                final_derived_scaled_noise_ratio = derived_scaled_noise_ratio
            else:
                final_derived_scaled_noise_ratio = ma.vstack(
                    (final_derived_scaled_noise_ratio, derived_scaled_noise_ratio)
                )
            if final_derived_scaled_noise_ratio_cut is None:
                final_derived_scaled_noise_ratio_cut = derived_scaled_noise_ratio_cut
            else:
                final_derived_scaled_noise_ratio_cut = ma.vstack(
                    (
                        final_derived_scaled_noise_ratio_cut,
                        derived_scaled_noise_ratio_cut,
                    )
                )
            if final_masked_spec_array_to_correct is None:
                final_masked_spec_array_to_correct = masked_spec_array_to_correct
            else:
                final_masked_spec_array_to_correct = ma.vstack(
                    (final_masked_spec_array_to_correct, masked_spec_array_to_correct)
                )
        if corrected_specs is None:
            continue
            # corrected_specs = original_specs
        pyclass.comm("sic message class s+i")

        if config["dump_original_spec_array"]:
            if not os.path.exists(".dump_arrays"):
                subprocess.call(["mkdir", ".dump_arrays"])
            file_name = ".dump_arrays/array_orignal_spectra_{}_{}_{}_{}.txt".format(
                mission_id, telescope, scan, subscan
            )
            np.savetxt(file_name, np.asarray(original_specs))

        if config["output_folder"]:
            image_path = "{}/{}_{}_{}_{}_{}_summary.png".format(
                config["output_folder"],
                mission_id,
                telescope,
                scan,
                subscan,
                config["pca_source"],
            )
            title = "{}_{}_{}_{}_{}".format(
                telescope, scan, subscan, config["pca_source"], config["mission_id"]
            )
            suptitle = (
                "mission_id: {} -- telescope: {} -- "
                "scan: {} -- subscan: {} -- pca_source: {}"
            ).format(mission_id, telescope, scan, subscan, config["pca_source"])
            pyclass.comm("set unit v")
            try:
                print(pyclass.gdict.rx.shape, pyclass.gdict.rx.__sicdata__)
                plot_summary(
                    original_spectra=original_specs,
                    corrected_spectra=corrected_specs,
                    pca_comp=decomposition.pca_comp,
                    sky_diff_array=decomposition.input_spectra.spectra,
                    pca=decomposition.pca,
                    image_path=image_path,
                    title=title,
                    suptitle=suptitle,
                    number_components=decomposition.number_components_used,
                    x_axis=pyclass.gdict.rx.__sicdata__,
                    cutoff=config["cutoff"],
                    channels_exclude=channels_exclude,
                    bad_value=bad_value,
                    good_channels_orig=good_channels,
                    bad_channels_orig=bad_channels,
                    good_channels=good_channel_array,
                    bad_channels=bad_channel_array,
                    cut_coefficients=config["cut_coefficients"],
                    derived_component_coefficients=derived_component_coefficients,
                    derived_scaled_noise_ratio=final_derived_scaled_noise_ratio,
                    derived_scaled_noise_ratio_cut=final_derived_scaled_noise_ratio_cut,
                    example_correction=example_correction,
                    global_structure=decomposition.global_structure,
                    global_derived_scaled_noise_ratio=global_derived_scaled_noise_ratio,
                    contoured_plot=final_contoured_plot,
                    masked_plot=final_masked_spec_array_to_correct,
                    first_kde_gauss_peak_params=first_kde_gauss_peak_params,
                    config=config,
                )
            except TypeError as e:
                print(e)
                raise SystemExit

        del original_specs
        del corrected_specs
        del masked_spec_array_to_correct
        # del contoured_plot
        del example_correction
        del spectrum_collection
        pyclass.comm("find")
        gc.collect()

    if config["export_components"]:
        telescope_groups = df.groupby(["telescope", "mission_id"])

        for name, group in telescope_groups:
            telescope = name[0]
            mission_id = name[1]
            config["telescope"] = telescope
            config["mission_id"] = mission_id
            pyclass.message(pyclass.seve.i, "PCA", "Writing components to outputfile ")
            # pca_comp_output_file = create_pickle_file_name(
            #     config
            #     )
            for n, spec in enumerate(decomposition.pca.components):
                pyclass.gdict.ry = spec
                pyclass.comm("modify source pca_comp_{}".format(n))
                pyclass.comm("modify telescope {}".format(telescope))
                pyclass.comm("write")

    if config["component_map"]:
        df = pd.DataFrame(component_map)
        df.to_csv(config["component_map"])
    if config["component_map"]:
        df = pd.DataFrame(component_map)
        df.to_csv(config["component_map"])
