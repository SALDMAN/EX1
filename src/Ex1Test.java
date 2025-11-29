import org.junit.jupiter.api.Test;
import java.util.Random;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as defined in Ex1.
 * @author Yair
 */
class Ex1Test {
    static final double[] P1 = {2, 0, 3, -1, 0}, P2 = {0.1, 0, 1, 0.1, 3};
    static double[] po1 = {2, 2}, po2 = {-3, 0.61, 0.2};
    static double[] po3 = {2, 1, -0.7, -0.02, 0.02};
    static double[] po4 = {-3, 0.61, 0.2};

    /** Tests the f(x) function for simple evaluation at x=0,1,2 */
    @Test
    void testF() {
        double fx0 = Ex1.f(po1, 0);
        double fx1 = Ex1.f(po1, 1);
        double fx2 = Ex1.f(po1, 2);
        assertEquals(fx0, 2, Ex1.EPS);
        assertEquals(fx1, 4, Ex1.EPS);
        assertEquals(fx2, 6, Ex1.EPS);
    }

    /** Tests that f(x) of sum of polynomials equals sum of f(x) values */
    @Test
    void testF2() {
        double x = Math.PI;
        double[] po12 = Ex1.add(po1, po2);
        double f1x = Ex1.f(po1, x);
        double f2x = Ex1.f(po2, x);
        double f12x = Ex1.f(po12, x);
        assertEquals(f1x + f2x, f12x, Ex1.EPS);
    }

    /** Tests addition and subtraction to recover original polynomial */
    @Test
    void testAdd() {
        double[] p12 = Ex1.add(po1, po2);
        double[] minus1 = {-1};
        double[] pp2 = Ex1.mul(po2, minus1);
        double[] p1 = Ex1.add(p12, pp2);
        assertTrue(Ex1.equals(p1, po1));
    }

    /** Tests that addition is commutative */
    @Test
    void testAdd2() {
        double[] p12 = Ex1.add(po1, po2);
        double[] p21 = Ex1.add(po2, po1);
        assertTrue(Ex1.equals(p12, p21));
    }

    /** Tests that adding ZERO polynomial returns the same polynomial */
    @Test
    void testAdd3() {
        double[] p1 = Ex1.add(po1, Ex1.ZERO);
        assertTrue(Ex1.equals(p1, po1));
    }

    /** Tests addition of polynomials of different lengths produces correct length */
    @Test
    void testAdd4() {
        double[] p1 = {2, 4, 3, 2};
        double[] p2 = {3, 2};
        double[] p3 = Ex1.add(p1, p2);
        if (p3.length != 4) {
            fail();
        }
    }

    /** Tests addition result equals expected coefficients */
    @Test
    void testAdd5() {
        double[] p1 = {2, 4, 3, 2};
        double[] p2 = {3, 2};
        double[] p3 = Ex1.add(p1, p2);
        double[] expected = {5, 6, 3, 2};
        assertTrue(Ex1.equals(p3, expected));
    }

    /** Tests multiplication of a polynomial by ZERO returns ZERO */
    @Test
    void testMul1() {
        double[] p1 = Ex1.mul(po1, Ex1.ZERO);
        assertTrue(Ex1.equals(p1, Ex1.ZERO));
    }

    /** Tests multiplication is commutative */
    @Test
    void testMul2() {
        double[] p12 = Ex1.mul(po1, po2);
        double[] p21 = Ex1.mul(po2, po1);
        assertTrue(Ex1.equals(p12, p21));
    }

    /** Tests that multiplication matches pointwise evaluation */
    @Test
    void testMulDoubleArrayDoubleArray() {
        double[] xx = {0, 1, 2, 3, 4.1, -15.2222};
        double[] p12 = Ex1.mul(po1, po2);
        for (int i = 0; i < xx.length; i++) {
            double x = xx[i];
            double f1x = Ex1.f(po1, x);
            double f2x = Ex1.f(po2, x);
            double f12x = Ex1.f(p12, x);
            assertEquals(f12x, f1x * f2x, Ex1.EPS);
        }
    }

    /** Tests derivative of polynomial reduces degree correctly */
    @Test
    void testDerivativeArrayDoubleArray() {
        double[] p = {1, 2, 3};
        double[] pt = {2, 6};
        double[] dp1 = Ex1.derivative(p);
        double[] dp2 = Ex1.derivative(dp1);
        double[] dp3 = Ex1.derivative(dp2);
        double[] dp4 = Ex1.derivative(dp3);
        assertTrue(Ex1.equals(dp1, pt));
        assertTrue(Ex1.equals(Ex1.ZERO, dp3));
        assertTrue(Ex1.equals(dp4, dp3));
    }

    /** Tests converting polynomial to string and back from string */
    @Test
    public void testFromString() {
        double[] p = {-1.1, 2.3, 3.1};
        String sp2 = "3.1x^2 +2.3x -1.1";
        String sp = Ex1.poly(p);
        double[] p1 = Ex1.getPolynomFromString(sp);
        double[] p2 = Ex1.getPolynomFromString(sp2);
        boolean isSame1 = Ex1.equals(p1, p);
        boolean isSame2 = Ex1.equals(p2, p);
        if (!isSame1) fail();
        if (!isSame2) fail();
        assertEquals(sp, Ex1.poly(p1));
    }

    /** Tests equals method for polynomials with small differences */
    @Test
    public void testEquals() {
        double[][] d1 = {{0}, {1}, {1, 2, 0, 0}};
        double[][] d2 = {Ex1.ZERO, {1 + Ex1.EPS / 2}, {1, 2}};
        double[][] xx = {{-2 * Ex1.EPS}, {1 + Ex1.EPS * 1.2}, {1, 2, Ex1.EPS / 2}};
        for (int i = 0; i < d1.length; i++) {
            assertTrue(Ex1.equals(d1[i], d2[i]));
        }
        for (int i = 0; i < d1.length; i++) {
            assertFalse(Ex1.equals(d1[i], xx[i]));
        }
    }

    /** Tests that sameValue returns consistent crossing point between polynomials */
    @Test
    public void testSameValue2() {
        double x1 = -4, x2 = 0;
        double rs1 = Ex1.sameValue(po1, po2, x1, x2, Ex1.EPS);
        double rs2 = Ex1.sameValue(po2, po1, x1, x2, Ex1.EPS);
        assertEquals(rs1, rs2, Ex1.EPS);
    }

    /** Tests area computation between two polynomials is symmetric */
    @Test
    public void testArea() {
        double x1 = -4, x2 = 0;
        double a1 = Ex1.area(po1, po2, x1, x2, 100);
        double a2 = Ex1.area(po2, po1, x1, x2, 100);
        assertEquals(a1, a2, Ex1.EPS);
    }

    /** Tests area calculation with different number of segments converges */
    @Test
    public void testArea2() {
        double[] po_a = Ex1.ZERO;
        double[] po_b = {0, 1};
        double x1 = -1;
        double x2 = 2;
        double a1 = Ex1.area(po_a, po_b, x1, x2, 1);
        double a2 = Ex1.area(po_a, po_b, x1, x2, 2);
        double a3 = Ex1.area(po_a, po_b, x1, x2, 3);
        double a100 = Ex1.area(po_a, po_b, x1, x2, 100);
        double area = 2.5;
        assertEquals(a1, area, Ex1.EPS);
        assertEquals(a2, area, Ex1.EPS);
        assertEquals(a3, area, Ex1.EPS);
        assertEquals(a100, area, Ex1.EPS);
    }

    /** Tests area calculation for two polynomials with known intersection */
    @Test
    public void testArea3() {
        double[] po_a = {2, 1, -0.7, -0.02, 0.02};
        double[] po_b = {6, 0.1, -0.2};
        double x1 = Ex1.sameValue(po_a, po_b, -10, -5, Ex1.EPS);
        double a1 = Ex1.area(po_a, po_b, x1, 6, 8);
        double area = 58.5658;
        assertEquals(a1, area, Ex1.EPS);
    }

    /** Tests derivative with specific coefficients */
    @Test
    public void testDerivative() {
        double[] p1 = {3, 2, 6, 4};
        double[] p2 = {2, 12, 12};
        double[] ans = Ex1.derivative(p1);
        for (int i = 0; i < ans.length; i++) {
            if (ans[i] != p2[i]) fail();
        }
    }

    /** Tests poly() string formatting of polynomial */
    @Test
    void testPoly() {
        double[] p1 = {7.0, -1.2, 0.0, -1.0, 3.5};
        String expected = "3.50x^4 - x^3 - 1.20x + 7.00";
        assertEquals(expected, Ex1.poly(p1), "Test Case 1 failed: Standard polynomial format error.");
    }

    /** Tests poly() string formatting of more complex polynomial */
    @Test
    void testPoly2() {
        double[] p3 = {3.33, -0.99, 5.55, -1.0, 0.0, 1.0, -2.55};
        String expected = "-2.55x^6 + x^5 - x^3 + 5.55x^2 - 0.99x + 3.33";
        assertEquals(expected, Ex1.poly(p3), "Test Case 3 failed: Complicated polynomial format error.");
    }

    /** Tests length calculation of a linear polynomial */
    @Test
    void testLengthLinear() {
        double[] p = {0, 2}; // y = 2x
        double len = Ex1.length(p, 0, 3, 3);
        double dx = (3.0 - 0.0) / 3;
        double dy = 2 * dx;
        double expected = Math.sqrt(dx*dx + dy*dy) * 3;
        assertEquals(expected, len, Ex1.EPS, "Length of linear function incorrect");
    }

    /** Tests length calculation of a quadratic polynomial */
    @Test
    void testLengthQuadratic() {
        double[] p = {0, 0, 1}; // y = x^2
        double len = Ex1.length(p, 0, 1, 1000);
        assertEquals(1.4789, len, 0.01, "Length of quadratic function incorrect");
    }

    /** Tests creation of linear polynomial from two points */
    @Test
    void testPolynomFromPointsLinear() {
        double[] xx = {1, 3};
        double[] yy = {2, 6};
        double[] coeffs = Ex1.PolynomFromPoints(xx, yy);
        assertArrayEquals(new double[]{0, 2}, coeffs, 1e-6, "Linear polynomial from points incorrect");
    }

    /** Tests creation of quadratic polynomial from three points */
    @Test
    void testPolynomFromPointsQuadratic() {
        double[] xx = {0, 1, 2};
        double[] yy = {1, 6, 17};
        double[] coeffs = Ex1.PolynomFromPoints(xx, yy);
        assertArrayEquals(new double[]{1, 2, 3}, coeffs, 1e-6, "Quadratic polynomial from points incorrect");
    }

    /** Tests creation of quadratic polynomial with different points */
    @Test
    void testPolynomFromPointsQuadratic2() {
        double[] xx = {0, 1, 2};
        double[] yy = {6, 1, 17};
        double[] coeffs = Ex1.PolynomFromPoints(xx, yy);
        assertArrayEquals(new double[]{6,-15.5, 10.5}, coeffs, 1e-6, "Quadratic polynomial from points incorrect (test 2)");
    }

    /** Tests invalid points (same x-values) return null */
    @Test
    void testPolynomFromPointsInvalid() {
        double[] xx = {1, 1};
        double[] yy = {2, 2};
        assertNull(Ex1.PolynomFromPoints(xx, yy), "Should return null for invalid points (same x values)");
    }
}
