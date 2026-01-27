import unittest
import numpy as np
import pylab
from pyteomics import pylab_aux, mgf, mzml
from pyteomics.auxiliary import PyteomicsError

"""
Tests for pylab_aux module
"""



class PlotLineTest(unittest.TestCase):
    """Tests for plot_line function"""

    def setUp(self):
        pylab.clf()

    def tearDown(self):
        pylab.close('all')

    def test_plot_line_with_xlim(self):
        """Test plotting a line with explicit xlim"""
        result = pylab_aux.plot_line(2, 3, xlim=(0, 10))
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 1)

    def test_plot_line_without_xlim(self):
        """Test plotting a line without xlim (uses current axis limits)"""
        pylab.xlim(0, 10)
        result = pylab_aux.plot_line(1, 5)
        self.assertIsNotNone(result)

    def test_plot_line_negative_slope(self):
        """Test plotting a line with negative slope"""
        result = pylab_aux.plot_line(-0.5, 10, xlim=(0, 20))
        self.assertIsNotNone(result)

    def test_plot_line_zero_intercept(self):
        """Test plotting a line through origin"""
        result = pylab_aux.plot_line(2, 0, xlim=(0, 5))
        self.assertIsNotNone(result)


class ScatterTrendTest(unittest.TestCase):
    """Tests for scatter_trend function"""

    def setUp(self):
        pylab.clf()
        self.x = np.array([1, 2, 3, 4, 5])
        self.y = np.array([2, 4, 5, 4, 5])

    def tearDown(self):
        pylab.close('all')

    def test_scatter_trend_with_separate_arrays(self):
        """Test scatter_trend with separate x and y arrays"""
        result = pylab_aux.scatter_trend(self.x, self.y)
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 4)  # (scatter, line, sigmas, legend)

    def test_scatter_trend_with_2d_array(self):
        """Test scatter_trend with 2D array"""
        data = np.column_stack([self.x, self.y])
        result = pylab_aux.scatter_trend(data)
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 4)

    def test_scatter_trend_no_trend_line(self):
        """Test scatter_trend without plotting trend line"""
        result = pylab_aux.scatter_trend(self.x, self.y, plot_trend=False)
        self.assertIsNone(result[1])  # line should be None

    def test_scatter_trend_with_sigmas(self):
        """Test scatter_trend with sigma lines"""
        result = pylab_aux.scatter_trend(self.x, self.y, plot_sigmas=True)
        self.assertIsNotNone(result[2])  # sigma lines should not be None

    def test_scatter_trend_no_legend(self):
        """Test scatter_trend without legend"""
        result = pylab_aux.scatter_trend(self.x, self.y, show_legend=False)
        self.assertIsNone(result[3])  # legend should be None

    def test_scatter_trend_with_title_and_labels(self):
        """Test scatter_trend with title and axis labels"""
        result = pylab_aux.scatter_trend(
            self.x, self.y,
            title='Test',
            xlabel='X axis',
            ylabel='Y axis'
        )
        self.assertIsNotNone(result)

    def test_scatter_trend_custom_regression(self):
        """Test scatter_trend with custom regression function"""
        def dummy_regression(x, y):
            return 1, 0, 0.9, 0.1
        result = pylab_aux.scatter_trend(self.x, self.y, regression=dummy_regression)
        self.assertIsNotNone(result)


class PlotFunction3DTest(unittest.TestCase):
    """Tests for plot_function_3d function"""

    def setUp(self):
        pylab.clf()
        self.x = np.linspace(0, 10, 5)
        self.y = np.linspace(0, 10, 5)
        self.function = lambda x, y: x + y

    def tearDown(self):
        pylab.close('all')

    def test_plot_function_3d_surface(self):
        """Test 3D surface plot"""
        pylab_aux.plot_function_3d(self.x, self.y, self.function, plot_type='surface')
        # Just verify it doesn't raise an exception

    def test_plot_function_3d_wireframe(self):
        """Test 3D wireframe plot"""
        pylab_aux.plot_function_3d(self.x, self.y, self.function, plot_type='wireframe')

    def test_plot_function_3d_scatter(self):
        """Test 3D scatter plot"""
        pylab_aux.plot_function_3d(self.x, self.y, self.function, plot_type='scatter')

    def test_plot_function_3d_contour(self):
        """Test 3D contour plot"""
        pylab_aux.plot_function_3d(self.x, self.y, self.function, plot_type='contour')

    def test_plot_function_3d_contourf(self):
        """Test 3D filled contour plot"""
        pylab_aux.plot_function_3d(self.x, self.y, self.function, plot_type='contourf')

    def test_plot_function_3d_invalid_plot_type(self):
        """Test 3D plot with invalid plot type"""
        with self.assertRaises(PyteomicsError):
            pylab_aux.plot_function_3d(self.x, self.y, self.function, plot_type='invalid')

    def test_plot_function_3d_with_labels(self):
        """Test 3D plot with axis labels"""
        pylab_aux.plot_function_3d(
            self.x, self.y, self.function,
            xlabel='X', ylabel='Y', zlabel='Z', title='3D Function'
        )


class PlotFunctionContourTest(unittest.TestCase):
    """Tests for plot_function_contour function"""

    def setUp(self):
        pylab.clf()
        self.x = np.linspace(0, 10, 5)
        self.y = np.linspace(0, 10, 5)
        self.function = lambda x, y: x**2 + y**2

    def tearDown(self):
        pylab.close('all')

    def test_plot_function_contour_filled(self):
        """Test filled contour plot"""
        pylab_aux.plot_function_contour(self.x, self.y, self.function, filling=True)

    def test_plot_function_contour_not_filled(self):
        """Test unfilled contour plot"""
        pylab_aux.plot_function_contour(self.x, self.y, self.function, filling=False)

    def test_plot_function_contour_custom_levels(self):
        """Test contour plot with custom number of levels"""
        pylab_aux.plot_function_contour(self.x, self.y, self.function, num_contours=20)

    def test_plot_function_contour_with_labels(self):
        """Test contour plot with labels"""
        pylab_aux.plot_function_contour(
            self.x, self.y, self.function,
            xlabel='X', ylabel='Y', title='Contour'
        )


class PlotQvalueCurveTest(unittest.TestCase):
    """Tests for plot_qvalue_curve function"""

    def setUp(self):
        pylab.clf()
        self.qvalues = np.array([0.001, 0.002, 0.005, 0.01, 0.02, 0.05])

    def tearDown(self):
        pylab.close('all')

    def test_plot_qvalue_curve_basic(self):
        """Test basic q-value curve plot"""
        result = pylab_aux.plot_qvalue_curve(self.qvalues)
        self.assertIsNotNone(result)

    def test_plot_qvalue_curve_custom_labels(self):
        """Test q-value curve with custom labels"""
        result = pylab_aux.plot_qvalue_curve(
            self.qvalues,
            xlabel='Custom q-value',
            ylabel='Custom PSMs',
            title='Custom Title'
        )
        self.assertIsNotNone(result)

    def test_plot_qvalue_curve_with_style(self):
        """Test q-value curve with custom style"""
        result = pylab_aux.plot_qvalue_curve(self.qvalues, 'r--', linewidth=2)
        self.assertIsNotNone(result)


class PlotSpectrumTest(unittest.TestCase):
    """Tests for plot_spectrum function"""

    def setUp(self):
        pylab.clf()
        self.spectrum = {
            'm/z array': np.array([100, 150, 200, 250, 300]),
            'intensity array': np.array([10, 50, 100, 30, 20]),
            'params': {
                'charge': [2],
                'pepmass': (250.4, 1100)
            }
        }

    def tearDown(self):
        pylab.close('all')

    def test_plot_spectrum_default_backend_centroided(self):
        """Test spectrum plot with default backend, centroided"""
        result = pylab_aux.plot_spectrum(self.spectrum, backend='default', centroided=True)
        self.assertIsNotNone(result)

    def test_plot_spectrum_default_backend_profile(self):
        """Test spectrum plot with default backend, profile"""
        result = pylab_aux.plot_spectrum(self.spectrum, backend='default', centroided=False)
        self.assertIsNotNone(result)

    def test_plot_spectrum(self):
        for backend in ['default', 'spectrum_utils']:
            with self.subTest(backend=backend):
                result = pylab_aux.plot_spectrum(self.spectrum, backend=backend)
                self.assertIsNotNone(result)

    def test_plot_spectrum_custom_labels(self):
        """Test spectrum plot with custom labels"""
        result = pylab_aux.plot_spectrum(
            self.spectrum,
            xlabel='m/z (Da)',
            ylabel='Intensity (counts)',
            title='Test Spectrum'
        )
        self.assertIsNotNone(result)

    def test_plot_spectrum_invalid_backend(self):
        """Test spectrum plot with invalid backend"""
        with self.assertRaises(PyteomicsError):
            pylab_aux.plot_spectrum(self.spectrum, backend='invalid_backend')


class AnnotateSpectrumTest(unittest.TestCase):
    """Tests for annotate_spectrum function"""

    def setUp(self):
        pylab.clf()
        self.spectrum = {
            'm/z array': np.array([100, 150, 200, 250, 300, 350]),
            'intensity array': np.array([10, 50, 100, 30, 20, 15]),
            'params': {
                'charge': [2],
                'pepmass': (250.5, 1000)
            }
        }
        self.peptide = 'PEPTIDE'

    def try_annotate_spectrum(self, **kwargs):
        for backend in ['default', 'spectrum_utils']:
            with self.subTest(backend=backend):
                kwargs['backend'] = backend
                result = pylab_aux.annotate_spectrum(
                    self.spectrum,
                    self.peptide,
                    **kwargs
                )
                self.assertIsNotNone(result)

    def tearDown(self):
        pylab.close('all')

    def test_annotate_spectrum(self):
        """Test spectrum annotation with default settings"""
        self.try_annotate_spectrum()

    def test_annotate_spectrum_custom_ion_types(self):
        """Test spectrum annotation with custom ion types"""
        ion_types = ['a', 'b', 'y']
        self.try_annotate_spectrum(ion_types=ion_types)

    def test_annotate_spectrum_custom_tolerance(self):
        """Test spectrum annotation with custom tolerance"""
        self.try_annotate_spectrum(ftol=0.5)

    def test_annotate_spectrum_override_charge(self):
        """Test spectrum annotation with overridden precursor charge"""
        self.try_annotate_spectrum(precursor_charge=3)

    def test_annotate_spectrum_custom_colors(self):
        """Test spectrum annotation with custom colors"""
        colors = {'b': '#FF0000', 'y': '#0000FF'}
        self.try_annotate_spectrum(colors=colors)

    def test_annotate_spectrum_with_title(self):
        """Test spectrum annotation with title"""
        self.try_annotate_spectrum(title='Annotated Spectrum')


class MirrorTest(unittest.TestCase):
    """Tests for mirror function"""

    def setUp(self):
        pylab.clf()
        self.spec_top = {
            'm/z array': np.array([100, 150, 200, 250, 300]),
            'intensity array': np.array([10, 50, 100, 30, 20]),
            'params': {
                'charge': [2],
                'pepmass': (250.4, 1100)}
        }
        self.spec_bottom = {
            'm/z array': np.array([100, 120, 200, 250, 300]),
            'intensity array': np.array([5, 40, 90, 25, 15]),
            'params': {
                'charge': [2],
                'pepmass': (250.5, 1000)}
        }
        self.peptide = 'PEPTIDE'

    def tearDown(self):
        pylab.close('all')

    def test_mirror_no_peptide_spectrum_utils_missing(self):
        """Test mirror plot without peptide (no annotation)"""
        if pylab_aux.sus is None:
            self.skipTest("spectrum_utils not available")
        result = pylab_aux.mirror(self.spec_top, self.spec_bottom, precursor_mz=250.5, precursor_charge=2)
        self.assertIsNotNone(result)

    def test_mirror_with_peptide_spectrum_utils_missing(self):
        """Test mirror plot with peptide (with annotation)"""
        if pylab_aux.sus is None:
            self.skipTest("spectrum_utils not available")
        result = pylab_aux.mirror(
            self.spec_top, self.spec_bottom,
            peptide=self.peptide,
            precursor_charge=2
        )
        self.assertIsNotNone(result)


class ExtractPrecursorInfoTest(unittest.TestCase):
    mgf_file = 'test.mgf'
    mzml_file = 'tiny.pwiz.1.1.mzML'

    def test_mgf_extraction(self):
        with mgf.IndexedMGF(self.mgf_file) as reader:
            spectrum = reader[0]
            charge = pylab_aux._get_precursor_charge(spectrum)
            mz = pylab_aux._get_precursor_mz(spectrum)
            self.assertEqual(charge, 2)
            self.assertAlmostEqual(mz, 983.6)

    def test_mzml_extraction(self):
        with mzml.MzML(self.mzml_file) as reader:
            spectrum = reader[1]
            charge = pylab_aux._get_precursor_charge(spectrum)
            mz = pylab_aux._get_precursor_mz(spectrum)
            self.assertEqual(charge, 2)
            self.assertAlmostEqual(mz, 445.34)

    def test_kwarg_override(self):
        spectrum = {}
        charge = pylab_aux._get_precursor_charge(spectrum, precursor_charge=3)
        mz = pylab_aux._get_precursor_mz(spectrum, precursor_mz=500.2)
        self.assertEqual(charge, 3)
        self.assertAlmostEqual(mz, 500.2)


if __name__ == '__main__':
    unittest.main()