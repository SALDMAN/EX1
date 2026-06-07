import java.io.*;

public class Kefel {

    public static void main(String[] argv) throws IOException {
        if (argv.length < 1) {
            System.err.println("Usage: java Kefel <k>");
            System.exit(1);
        }

        int k = Integer.parseInt(argv[0]);

        PrintWriter pw = new PrintWriter(new FileWriter("kefel.s"));

        pw.println(".section .text");
        pw.println(".globl kefel");
        pw.println("kefel:");

        writeMultiply(pw, k);

        pw.println("\tret");
        pw.close();
    }

    static void writeMultiply(PrintWriter pw, int k) {
        // k=1: just move
        if (k == 1) {
            pw.println("\tmovq %rdi,%rax");
            return;
        }

        // Rule 1: k = 3, 5, or 9 → use lea
        if (k == 3) {
            pw.println("\tleaq (%rdi,%rdi,2),%rax");
            return;
        }
        if (k == 5) {
            pw.println("\tleaq (%rdi,%rdi,4),%rax");
            return;
        }
        if (k == 9) {
            pw.println("\tleaq (%rdi,%rdi,8),%rax");
            return;
        }

        // Rule 2: power of 2 (single 1-bit) → single shift
        if (k > 0 && (k & (k - 1)) == 0) {
            int shift = Integer.numberOfTrailingZeros(k);
            pw.println("\tmovq %rdi,%rax");
            pw.println("\tshlq $" + shift + ",%rax");
            return;
        }

        // Check if k is a single run of consecutive 1-bits
        int highBit = 31 - Integer.numberOfLeadingZeros(k);
        int lowBit  = Integer.numberOfTrailingZeros(k);
        int runLength = highBit - lowBit + 1;
        int mask = ((1 << runLength) - 1) << lowBit;

        if (k == mask) {
            if (runLength == 2) {
                // Rule 3: 2 consecutive bits → add of two shifts
                pw.println("\tmovq %rdi,%rax");
                pw.println("\tshlq $" + (lowBit + 1) + ",%rax");
                if (lowBit == 0) {
                    pw.println("\taddq %rdi,%rax");
                } else {
                    pw.println("\tmovq %rdi,%rcx");
                    pw.println("\tshlq $" + lowBit + ",%rcx");
                    pw.println("\taddq %rcx,%rax");
                }
            } else {
                // Rule 4: 3+ consecutive bits → subtract of two shifts
                pw.println("\tmovq %rdi,%rax");
                pw.println("\tshlq $" + (highBit + 1) + ",%rax");
                if (lowBit == 0) {
                    pw.println("\tsubq %rdi,%rax");
                } else {
                    pw.println("\tmovq %rdi,%rcx");
                    pw.println("\tshlq $" + lowBit + ",%rcx");
                    pw.println("\tsubq %rcx,%rax");
                }
            }
            return;
        }

        // General case: multiple runs → sum individual shifts for each set bit
        boolean first = true;
        for (int i = 31; i >= 0; i--) {
            if ((k & (1 << i)) != 0) {
                if (first) {
                    pw.println("\tmovq %rdi,%rax");
                    if (i > 0) pw.println("\tshlq $" + i + ",%rax");
                    first = false;
                } else {
                    if (i == 0) {
                        pw.println("\taddq %rdi,%rax");
                    } else {
                        pw.println("\tmovq %rdi,%rcx");
                        pw.println("\tshlq $" + i + ",%rcx");
                        pw.println("\taddq %rcx,%rax");
                    }
                }
            }
        }
    }
}