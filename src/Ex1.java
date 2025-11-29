import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.lang.Math;

public class Ex1 {
    /** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
    public static final double EPS = 0.001;
    /** The zero polynomial function is represented as an array with a single (0) entry. */
    public static final double[] ZERO = {0};
    public enum AreaMode { SIGNED, SIGNED_ABS, INTEGRAL_ABS }


    public static double f(double[] poly, double x) {
        double ans = 0;
        for(int i=0;i<poly.length;i++) {
            double c = Math.pow(x, i);
            ans += c*poly[i];
        }
        return ans;
    }

    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p,x1);

        if (x2 - x1 < eps) {
            return (x1 + x2) / 2.0;
        }

        double x12 = (x1+x2)/2;
        double f12 = f(p,x12);

        if (Math.abs(f12) < eps) {
            return x12;
        }

        if (f1 * f12 <= 0) {
            return root_rec(p, x1, x12, eps);
        } else {
            return root_rec(p, x12, x2, eps);
        }
    }

    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        if (xx == null || yy == null || xx.length != yy.length || xx.length < 2 || xx.length > 3) {
            return null;
        }

        int n = xx.length;

        if (n == 2) {
            // Linear: y = a*x + b
            double x0 = xx[0], y0 = yy[0];
            double x1 = xx[1], y1 = yy[1];

            if (Math.abs(x0 - x1) < Ex1.EPS) return null;

            double a = (y1 - y0) / (x1 - x0);
            double b = y0 - a * x0;

            return new double[]{b, a}; // {constant, linear}
        }

        if (n == 3) {
            // Quadratic: y = a*x^2 + b*x + c
            double x0 = xx[0], y0 = yy[0];
            double x1 = xx[1], y1 = yy[1];
            double x2 = xx[2], y2 = yy[2];

            // Solve system using determinants (Cramer's rule)
            double det = x0*x0*(x1 - x2) - x1*x1*(x0 - x2) + x2*x2*(x0 - x1);
            if (Math.abs(det) < Ex1.EPS) return null;

            double a = (y0*(x1 - x2) - y1*(x0 - x2) + y2*(x0 - x1)) / det;
            double b = (y0*(x2*x2 - x1*x1) - y1*(x2*x2 - x0*x0) + y2*(x1*x1 - x0*x0)) / det;
            double c = (y0*(x1*x2*(x1 - x2)) - y1*(x0*x2*(x0 - x2)) + y2*(x0*x1*(x0 - x1))) / det;

            return new double[]{c, b, a}; // {constant, linear, quadratic}
        }

        return null;
    }

    public static boolean equals(double[] p1, double[] p2) {
        int n = Math.max(p1.length, p2.length) - 1;
        for (int i = 0; i <= n + 1; i++) {
            double x = i;
            if (Math.abs(f(p1, x) - f(p2, x)) > EPS) return false;
        }
        return true;
    }

    /**
     * Computes a String representing the polynomial function.
     */
    /**
     * Computes a String representing the polynomial function.
     */
    /**
     * Computes a String representing the polynomial function with correct signs.
     */
    /**
     * Computes a String representing the polynomial function with correct signs.
     */
    public static String poly(double[] poly) {
        if (poly == null || poly.length == 0) return "0.0";

        StringBuilder sb = new StringBuilder();
        boolean firstTerm = true;

        // עובר על המקדמים מהדרגה הגבוהה לנמוכה
        for (int i = poly.length - 1; i >= 0; i--) {
            double coef = poly[i];

            // דלג על מקדמים אפסיים (תוך שימוש ב-EPS המוגדר במחלקה)
            if (Math.abs(coef) < EPS) continue;

            // --- טיפול בסימן וברווח ---
            if (!firstTerm) {
                // אם זה לא האיבר הראשון:
                // אם המקדם חיובי, הוסף " + ".
                // אם המקדם שלילי, הוסף " - ".
                sb.append(coef >= 0 ? " + " : " - ");
            } else {
                // אם זה האיבר הראשון:
                // אם המקדם שלילי, הוסף "-". אחרת, לא מוסיפים סימן +.
                if (coef < 0) sb.append("-");
                firstTerm = false;
            }

            // --- טיפול בערך המקדם ---
            double abs = Math.abs(coef);

            // דרגה 0 (קבוע)
            if (i == 0) {
                sb.append(String.format("%.2f", abs));
            }
            // דרגה 1 ומעלה (עם x)
            else {
                // אם המקדם המוחלט שונה מ-1 (בטווח EPS), מדפיסים אותו.
                if (Math.abs(abs - 1.0) > EPS) {
                    sb.append(String.format("%.2f", abs));
                }

                sb.append("x");

                // דרגה > 1
                if (i > 1) {
                    sb.append("^").append(i);
                }
            }
        }

        // אין צורך ב-trim() מכיוון שטיפלנו ברווחים באופן מדויק.
        return firstTerm ? "0.0" : sb.toString();
    }

    public static double sameValue(double[] p1, double[] p2,
                                   double x1, double x2, double eps) {

        double prevX = x1;
        double prevVal = f(p1, prevX) - f(p2, prevX);

        int samples = 500000;         // דגימה מאוד צפופה
        double step = (x2 - x1) / samples;

        for (int i = 1; i <= samples; i++) {
            double x = x1 + i * step;
            double currVal = f(p1, x) - f(p2, x);

            if (prevVal * currVal <= 0) {
                // מצאנו קטע עם שינוי סימן -> עושים חיפוש בינארי על הקטע הזה
                return rootBinary(p1, p2, prevX, x, eps);
            }

            prevX = x;
            prevVal = currVal;
        }

        return Double.NaN;
    }



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
                // יש חיתוך בתוך הטרפז
                double root = rootBinary(p1, p2, left, right, EPS);
                // חישוב שטח בנפרד לכל חלק
                double area1 = (Math.abs(f(p1, left) - f(p2, left)) + Math.abs(f(p1, root) - f(p2, root))) / 2.0 * (root - left);
                double area2 = (Math.abs(f(p1, root) - f(p2, root)) + Math.abs(f(p1, right) - f(p2, right))) / 2.0 * (right - root);
                totalArea += area1 + area2;
            } else {
                // טרפז רגיל
                totalArea += (Math.abs(fLeft) + Math.abs(fRight)) / 2.0 * dx;
            }
        }

        return totalArea;
    }



    // rootBinary מדויק
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




    public static double[] getPolynomFromString(String p) {
        if (p == null || p.isEmpty()) return new double[]{0};

        String cleaned = p.replaceAll("\\s+", "");
        if (cleaned.isEmpty()) return new double[]{0};

        // החלפת - ב-+-, אך רק אם הוא לא בתחילת המחרוזת
        if (cleaned.startsWith("-")) {
            cleaned = "0" + cleaned;
        }
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

                if (coefStr.isEmpty()) coeff = 1.0;
                else coeff = Double.parseDouble(coefStr);
            }
            else if (cleanMonom.contains("x")) {
                degree = 1;
                int indexX = cleanMonom.indexOf("x");
                String coefStr = cleanMonom.substring(0, indexX);

                if (coefStr.isEmpty()) coeff = 1.0;
                else coeff = Double.parseDouble(coefStr);
            }
            else {
                degree = 0;
                coeff = Double.parseDouble(cleanMonom);
            }

            if (isNegative) coeff = -coeff;

            if (degree > maxDegree) maxDegree = degree;
            parsedMonoms.add(new double[]{coeff, degree});
        }

        double[] ans = new double[maxDegree + 1];
        for (double[] mono : parsedMonoms) {
            ans[(int) mono[1]] += mono[0];
        }

        // צמצום אפסים מובילים
        int effectiveLen = ans.length;
        while (effectiveLen > 1 && Math.abs(ans[effectiveLen - 1]) < EPS) {
            effectiveLen--;
        }

        if (effectiveLen < ans.length) {
            double[] trimmedAns = new double[effectiveLen];
            System.arraycopy(ans, 0, trimmedAns, 0, effectiveLen);
            return trimmedAns;
        }

        return ans;
    }

    public static double[] add(double[] p1, double[] p2) {
        int maxLen = Math.max(p1.length, p2.length);
        double [] ans = new double[maxLen];

        for (int i = 0; i < maxLen; i++) {
            double c1 = (i < p1.length) ? p1[i] : 0;
            double c2 = (i < p2.length) ? p2[i] : 0;
            ans[i] = c1 + c2;
        }

        // צמצום אפסים מובילים
        int effectiveLen = ans.length;
        while (effectiveLen > 1 && Math.abs(ans[effectiveLen - 1]) < EPS) {
            effectiveLen--;
        }

        // אם הצמצום הפך את המערך לפולינום האפס (אורך 1, ערך 0)
        if (effectiveLen == 1 && Math.abs(ans[0]) < EPS) {
            return ZERO;
        }

        if (effectiveLen < ans.length) {
            double[] trimmedAns = new double[effectiveLen];
            System.arraycopy(ans, 0, trimmedAns, 0, effectiveLen);
            return trimmedAns;
        }

        return ans;
    }

    public static double[] mul(double[] p1, double[] p2) {
        if (p1 == null || p2 == null || p1.length == 0 || p2.length == 0) return ZERO;

        double [] ans = new double[p1.length + p2.length - 1];

        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                ans[i+j] += p1[i] * p2[j];
            }
        }

        // צמצום אפסים מובילים (כפל אפס בפולינום צריך להחזיר רק {0})
        int effectiveLen = ans.length;
        while (effectiveLen > 1 && Math.abs(ans[effectiveLen - 1]) < EPS) {
            effectiveLen--;
        }

        if (effectiveLen < ans.length) {
            double[] trimmedAns = new double[effectiveLen];
            System.arraycopy(ans, 0, trimmedAns, 0, effectiveLen);
            return trimmedAns;
        }

        return ans;
    }

    public static double[] derivative(double[] po) {
        if (po.length <= 1) return new double[]{0};

        double[] ans = new double[po.length - 1];

        for (int i = 1; i < po.length; i++) {
            ans[i - 1] = po[i] * i;
        }

        return ans;
    }
}