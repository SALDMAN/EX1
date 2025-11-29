import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.lang.Math;

/**
 * A utility class for performing common mathematical operations on real-valued
 * polynomials represented by coefficient arrays.
 * This class is the second assignment counts from 0 in course Intro to computer science in Ariel University.
 * @author Yair
 */
public class Ex1 {

    /**
     * Epsilon value used as tolerance threshold for floating-point comparison.
     * Two values are considered equal if their absolute difference is below EPS.
     */
    public static final double EPS = 0.001;

    /**
     * Represents the zero polynomial, p(x) = 0.
     * Implemented as a single-element array containing 0.
     */
    public static final double[] ZERO = {0};

    /**
     * Enumeration of supported computation modes for area calculation.
     */
    public enum AreaMode { SIGNED, SIGNED_ABS, INTEGRAL_ABS }

    /**
     * Computes the value of a polynomial p(x) for a given x.
     *
     * @param poly The polynomial coefficients array.
     * @param x The x-value where f(x) is evaluated.
     * @return The computed value f(x).
     */
    public static double f(double[] poly, double x) {
        double ans = 0;
        for(int i = 0; i < poly.length; i++) {
            ans += poly[i] * Math.pow(x, i);
        }
        return ans;
    }

    /**
     * Recursively finds a root of the polynomial p(x) within interval [x1, x2]
     * using a binary search (bisection) method.
     *
     * @param p   The polynomial.
     * @param x1  Left boundary of interval.
     * @param x2  Right boundary of interval.
     * @param eps Stopping threshold for interval size or function value.
     * @return An approximate root of p(x) in the interval.
     */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p, x1);

        if (x2 - x1 < eps) {
            return (x1 + x2) / 2.0;
        }

        double mid = (x1 + x2) / 2;
        double fMid = f(p, mid);

        if (Math.abs(fMid) < eps) return mid;

        if (f1 * fMid <= 0) return root_rec(p, x1, mid, eps);
        else return root_rec(p, mid, x2, eps);
    }

    /**
     * Constructs a polynomial of degree 1 or 2 that exactly fits 2 or 3 given points.
     *
     * @param xx x-coordinates of given points.
     * @param yy y-coordinates of given points.
     * @return Polynomial coefficients array or null if invalid input.
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        if (xx == null || yy == null || xx.length != yy.length || xx.length < 2 || xx.length > 3) {
            return null;
        }

        int n = xx.length;

        if (n == 2) {
            double x0 = xx[0], y0 = yy[0];
            double x1 = xx[1], y1 = yy[1];

            if (Math.abs(x0 - x1) < Ex1.EPS) return null;

            double a = (y1 - y0) / (x1 - x0);
            double b = y0 - a * x0;

            return new double[]{b, a};
        }

        if (n == 3) {
            double x0 = xx[0], y0 = yy[0];
            double x1 = xx[1], y1 = yy[1];
            double x2 = xx[2], y2 = yy[2];

            double det = x0*x0*(x1 - x2) - x1*x1*(x0 - x2) + x2*x2*(x0 - x1);
            if (Math.abs(det) < Ex1.EPS) return null;

            double a = (y0*(x1 - x2) - y1*(x0 - x2) + y2*(x0 - x1)) / det;
            double b = (y0*(x2*x2 - x1*x1) - y1*(x2*x2 - x0*x0) + y2*(x1*x1 - x0*x0)) / det;
            double c = (y0*(x1*x2*(x1 - x2)) - y1*(x0*x2*(x0 - x2)) + y2*(x0*x1*(x0 - x1))) / det;

            return new double[]{c, b, a};
        }

        return null;
    }

    /**
     * Checks equality of two polynomials by evaluating them at several sample points.
     *
     * @param p1 First polynomial.
     * @param p2 Second polynomial.
     * @return true if polynomials represent the same function within tolerance.
     */
    public static boolean equals(double[] p1, double[] p2) {
        int n = Math.max(p1.length, p2.length) - 1;
        for (int i = 0; i <= n + 1; i++) {
            double x = i;
            if (Math.abs(f(p1, x) - f(p2, x)) > EPS) return false;
        }
        return true;
    }

    /**
     * Converts a polynomial into a human-readable algebraic string.
     *
     * @param poly Polynomial coefficients.
     * @return String representation of the polynomial.
     */
    public static String poly(double[] poly) {
        if (poly == null || poly.length == 0) return "0.0";

        StringBuilder sb = new StringBuilder();
        boolean firstTerm = true;

        for (int i = poly.length - 1; i >= 0; i--) {
            double coef = poly[i];
            if (Math.abs(coef) < EPS) continue;

            if (!firstTerm) {
                sb.append(coef >= 0 ? " + " : " - ");
            } else {
                if (coef < 0) sb.append("-");
                firstTerm = false;
            }

            double abs = Math.abs(coef);
            if (i == 0) sb.append(String.format("%.2f", abs));
            else {
                if (Math.abs(abs - 1.0) > EPS) sb.append(String.format("%.2f", abs));
                sb.append("x");
                if (i > 1) sb.append("^").append(i);
            }
        }

        return firstTerm ? "0.0" : sb.toString();
    }

    /**
     * Finds an x in [x1, x2] such that p1(x) == p2(x).
     *
     * @param p1 First polynomial.
     * @param p2 Second polynomial.
     * @param x1 Interval start.
     * @param x2 Interval end.
     * @param eps Accuracy threshold.
     * @return Intersection x-value or NaN if none found.
     */
    public static double sameValue(double[] p1, double[] p2,
                                   double x1, double x2, double eps) {

        double prevX = x1;
        double prevVal = f(p1, prevX) - f(p2, prevX);

        int samples = 500000;
        double step = (x2 - x1) / samples;

        for (int i = 1; i <= samples; i++) {
            double x = x1 + i * step;
            double currVal = f(p1, x) - f(p2, x);

            if (prevVal * currVal <= 0) {
                return rootBinary(p1, p2, prevX, x, eps);
            }

            prevX = x;
            prevVal = currVal;
        }

        return Double.NaN;
    }

    /**
     * Estimates the arc length of the polynomial curve between x1 and x2.
     *
     * @param p The polynomial.
     * @param x1 Start of interval.
     * @param x2 End of interval.
     * @param numberOfSegments Number of linear subdivisions.
     * @return Approximate arc length.
     */
    public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        if (numberOfSegments <= 0) throw new IllegalArgumentException("numberOfSegments must be positive");

        double dx = (x2 - x1) / numberOfSegments;
        double length = 0.0;

        double prevX = x1;
        double prevY = f(p, prevX);

        for (int i = 1; i <= numberOfSegments; i++) {
            double currX = x1 + i * dx;
            double currY = f(p, currX);

            double segment = Math.sqrt(Math.pow(dx, 2) + Math.pow(currY - prevY, 2));
            length += segment;

            prevX = currX;
            prevY = currY;
        }

        return length;
    }

    /**
     * Computes numeric area between two polynomials p1 and p2 over [x1, x2].
     *
     * @param p1 First polynomial.
     * @param p2 Second polynomial.
     * @param x1 Start of interval.
     * @param x2 End of interval.
     * @param n  Number of trapezoids.
     * @return Approximate area (always non-negative).
     */
    public static double area(double[] p1, double[] p2, double x1, double x2, int n) {
        if (n <= 0) throw new IllegalArgumentException("n must be positive");
        if (x1 > x2) { double t = x1; x1 = x2; x2 = t; }

        double dx = (x2 - x1) / n;
        double totalArea = 0.0;

        for (int i = 0; i < n; i++) {
            double left = x1 + i * dx;
            double right = left + dx;

            double fLeft = f(p1, left) - f(p2, left);
            double fRight = f(p1, right) - f(p2, right);

            if (fLeft * fRight < 0) {
                double root = rootBinary(p1, p2, left, right, EPS);
                double area1 = (Math.abs(f(p1, left) - f(p2, left)) + Math.abs(f(p1, root) - f(p2, root))) / 2.0 * (root - left);
                double area2 = (Math.abs(f(p1, root) - f(p2, root)) + Math.abs(f(p1, right) - f(p2, right))) / 2.0 * (right - root);
                totalArea += area1 + area2;
            } else {
                totalArea += (Math.abs(fLeft) + Math.abs(fRight)) / 2.0 * dx;
            }
        }

        return totalArea;
    }

    private static double rootBinary(double[] p1, double[] p2, double left, double right, double eps) {
        while (right - left > eps) {
            double mid = (left + right) / 2.0;
            double fLeft = f(p1, left) - f(p2, left);
            double fMid = f(p1, mid) - f(p2, mid);

            if (fLeft * fMid <= 0) right = mid;
            else left = mid;
        }
        return (left + right) / 2.0;
    }

    /**
     * Parses a polynomial string into coefficient array.
     *
     * @param p String representation of a polynomial.
     * @return Coefficient array.
     */
    public static double[] getPolynomFromString(String p) {
        if (p == null || p.isEmpty()) return new double[]{0};

        String cleaned = p.replaceAll("\\s+", "");
        if (cleaned.isEmpty()) return new double[]{0};

        if (cleaned.startsWith("-")) cleaned = "0" + cleaned;
        cleaned = cleaned.replace("-", "+-");
        String[] monoms = cleaned.split("\\+");

        int maxDegree = 0;
        List<double[]> parsedMonoms = new ArrayList<>();

        for (String m : monoms) {
            if (m.isEmpty() || m.equals("0")) continue;

            int degree = 0;
            double coeff = 0.0;
            String cleanMonom = m;

            boolean isNegative = m.startsWith("-");
            if (isNegative) cleanMonom = m.substring(1);

            if (cleanMonom.contains("x^")) {
                int indexX = cleanMonom.indexOf("x");
                degree = Integer.parseInt(cleanMonom.substring(cleanMonom.indexOf("^") + 1));
                String coefStr = cleanMonom.substring(0, indexX);
                coeff = coefStr.isEmpty() ? 1.0 : Double.parseDouble(coefStr);
            } else if (cleanMonom.contains("x")) {
                degree = 1;
                String coefStr = cleanMonom.substring(0, cleanMonom.indexOf("x"));
                coeff = coefStr.isEmpty() ? 1.0 : Double.parseDouble(coefStr);
            } else {
                degree = 0;
                coeff = Double.parseDouble(cleanMonom);
            }

            if (isNegative) coeff = -coeff;

            if (degree > maxDegree) maxDegree = degree;
            parsedMonoms.add(new double[]{coeff, degree});
        }

        double[] ans = new double[maxDegree + 1];
        for (double[] mono : parsedMonoms) ans[(int) mono[1]] += mono[0];

        int effectiveLen = ans.length;
        while (effectiveLen > 1 && Math.abs(ans[effectiveLen - 1]) < EPS) effectiveLen--;

        if (effectiveLen < ans.length) {
            double[] trimmedAns = new double[effectiveLen];
            System.arraycopy(ans, 0, trimmedAns, 0, effectiveLen);
            return trimmedAns;
        }

        return ans;
    }

    /**
     * Adds two polynomials.
     *
     * @param p1 First polynomial.
     * @param p2 Second polynomial.
     * @return Sum polynomial.
     */
    public static double[] add(double[] p1, double[] p2) {
        int maxLen = Math.max(p1.length, p2.length);
        double [] ans = new double[maxLen];

        for (int i = 0; i < maxLen; i++) {
            double c1 = (i < p1.length) ? p1[i] : 0;
            double c2 = (i < p2.length) ? p2[i] : 0;
            ans[i] = c1 + c2;
        }

        int effectiveLen = ans.length;
        while (effectiveLen > 1 && Math.abs(ans[effectiveLen - 1]) < EPS) effectiveLen--;

        if (effectiveLen == 1 && Math.abs(ans[0]) < EPS) return ZERO;

        if (effectiveLen < ans.length) {
            double[] trimmedAns = new double[effectiveLen];
            System.arraycopy(ans, 0, trimmedAns, 0, effectiveLen);
            return trimmedAns;
        }

        return ans;
    }

    /**
     * Multiplies two polynomials.
     *
     * @param p1 First polynomial.
     * @param p2 Second polynomial.
     * @return Product polynomial.
     */
    public static double[] mul(double[] p1, double[] p2) {
        if (p1 == null || p2 == null || p1.length == 0 || p2.length == 0) return ZERO;

        double [] ans = new double[p1.length + p2.length - 1];

        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                ans[i+j] += p1[i] * p2[j];
            }
        }

        int effectiveLen = ans.length;
        while (effectiveLen > 1 && Math.abs(ans[effectiveLen - 1]) < EPS) effectiveLen--;

        if (effectiveLen < ans.length) {
            double[] trimmedAns = new double[effectiveLen];
            System.arraycopy(ans, 0, trimmedAns, 0, effectiveLen);
            return trimmedAns;
        }

        return ans;
    }

    /**
     * Computes the derivative polynomial of p(x).
     *
     * @param po Polynomial coefficients.
     * @return The derivative polynomial.
     */
    public static double[] derivative(double[] po) {
        if (po.length <= 1) return new double[]{0};

        double[] ans = new double[po.length - 1];

        for (int i = 1; i < po.length; i++) {
            ans[i - 1] = po[i] * i;
        }

        return ans;
    }
}
