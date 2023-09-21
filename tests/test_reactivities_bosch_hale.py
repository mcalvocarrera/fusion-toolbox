import pytest
from fusion_toolbox.reactivity import BoschHaleReactivity

@pytest.mark.parametrize(
    "test, expected, precision",
    [
        (1.0, 6.857e-21, 1e-21),
        (2.0, 2.977e-19, 1e-19),
        (5.0, 1.366e-17, 1e-18),
        (10.0, 1.136e-16, 1e-17),
        (20.0, 4.330e-16, 1e-16),
    ])
def test_DTn_reactivity(test, expected, precision):
    """Test of the DTn fusion reactivity
    Tested according to table VIII in H.-S. Bosch and G.M. Hale 1992 Nucl. Fusion 32 611

    Parameters
    ---------
    test : float
        Input ion temperature for the test
    expected : float
        Expected answer from table VIII in Bosch-Hale '92
    precision : int 
        Precision for the comparision between the Bosch-Hale parameterized value
        and the tabulated value
    """
    DTn_r = BoschHaleReactivity("DTn")
    assert DTn_r.get_reactivity(test) == pytest.approx(expected, abs=precision)

@pytest.mark.parametrize(
    "test, expected, precision",
    [
        (1.0, 9.933e-23, 1e-24),
        (2.0, 3.110e-21, 1e-22),
        (5.0, 9.128e-20, 1e-21),
        (10.0, 6.023e-19, 1e-20),
        (20.0, 2.603e-18, 1e-19),
    ])
def test_DDn_reactivity(test, expected, precision):
    """Test of the DDn fusion reactivity
    Tested according to table VIII in H.-S. Bosch and G.M. Hale 1992 Nucl. Fusion 32 611

    Parameters
    ---------
    test : float
        Input ion temperature for the test
    expected : float
        Expected answer from table VIII in Bosch-Hale '92
    precision : int 
        Precision for the comparision between the Bosch-Hale parameterized value
        and the tabulated value
    """
    DDn_r = BoschHaleReactivity("DDn")
    assert DDn_r.get_reactivity(test) == pytest.approx(expected, abs=precision)

@pytest.mark.parametrize(
    "test, expected, precision",
    [
        (1.0, 1.017e-22, 1e-23),
        (2.0, 3.150e-21, 1e-22),
        (5.0, 9.024e-20, 1e-21),
        (10.0, 5.781e-19, 1e-20),
        (20.0, 2.399e-18, 1e-20),
    ])
def test_DDp_reactivity(test, expected, precision):
    """Test of the DDp fusion reactivity
    Tested according to table VIII in H.-S. Bosch and G.M. Hale 1992 Nucl. Fusion 32 611

    Parameters
    ---------
    test : float
        Input ion temperature for the test
    expected : float
        Expected answer from table VIII in Bosch-Hale '92
    precision : int 
        Precision for the comparision between the Bosch-Hale parameterized value
        and the tabulated value
    """
    DDp_r = BoschHaleReactivity("DDp")
    assert DDp_r.get_reactivity(test) == pytest.approx(expected, abs=precision)

@pytest.mark.parametrize(
    "test, expected, precision",
    [
        (1.0, 3.057e-26, 1e-27),
        (2.0, 1.399e-23, 1e-24),
        (5.0, 6.377e-21, 1e-22),
        (10.0, 2.126e-19, 1e-20),
        (20.0, 3.482e-18, 1e-21),
    ])
def test_D3Hep_reactivity(test, expected, precision):
    """Test of the D3Hep fusion reactivity
    Tested according to table VIII in H.-S. Bosch and G.M. Hale 1992 Nucl. Fusion 32 611

    Parameters
    ---------
    test : float
        Input ion temperature for the test
    expected : float
        Expected answer from table VIII in Bosch-Hale '92
    precision : int 
        Precision for the comparision between the Bosch-Hale parameterized value
        and the tabulated value
    """
    D3Hep_r = BoschHaleReactivity("D3Hep")
    assert D3Hep_r.get_reactivity(test) == pytest.approx(expected, abs=precision)

def test_wrong_reaction_exception():
    with pytest.raises(ValueError):
        DTgamma_r = BoschHaleReactivity("DTgamma")
