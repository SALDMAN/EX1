import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Kefel {
    public static void main(String[] args) {
        if (args.length < 1) {
            System.out.println("Usage: java Kefel <constant>");
            return;
        }

        int k = Integer.parseInt(args[0]);

        try (PrintWriter out = new PrintWriter(new FileWriter("kefel.s"))) {
            // כתיבת הפתיח הקבוע של קובץ האסמבלי
            out.println(".section .text");
            out.println(".globl kefel");
            out.println("kefel:");

            // טיפול במקרה קצה שבו k שווה ל-0
            if (k == 0) {
                out.println("movq $0, %rax");
                out.println("ret");
                return;
            }

            // כלל 1: אם k הוא 3, 5, או 9 - שימוש בפקודת lea בודדת
            if (k == 3 || k == 5 || k == 9) {
                int scale = k - 1;
                out.println("leaq (%rdi,%rdi," + scale + "), %rax");
                out.println("ret");
                return;
            }

            // ניתוח המבנה הבינארי של המספר לצורך שאר הכללים
            int lowestSetBit = Integer.numberOfTrailingZeros(k);
            int highestSetBit = 31 - Integer.numberOfLeadingZeros(k);
            int bitCount = Integer.bitCount(k);
            int blockLength = highestSetBit - lowestSetBit + 1;

            // כלל 2: אם k מכיל ביט אחד בלבד של 1 (חזקה של 2) -> הזזה אחת
            if (bitCount == 1) {
                out.println("movq %rdi, %rax");
                if (lowestSetBit > 0) {
                    out.println("shlq $" + lowestSetBit +, %rax");
                }
                out.println("ret");
                return;
            }

            // כלל 3: אם k מכיל בדיוק 2 ביטים של 1 (לא משנה אם רצופים או לא, להשגת מינימום שורות) -> חיבור של שתי הזזות
            if (bitCount == 2) {
                int firstBit = lowestSetBit;
                int secondBit = 31 - Integer.numberOfLeadingZeros(k ^ (1 << firstBit));

                out.println("movq %rdi, %rax");
                if (secondBit > 0) {
                    out.println("shlq $" + secondBit + ", %rax");
                }
                out.println("movq %rdi, %rcx");
                if (firstBit > 0) {
                    out.println("shlq $" + firstBit + ", %rcx");
                }
                out.println("addq %rcx, %rax");
                out.println("ret");
                return;
            }

            // כלל 4: אם k הוא רצף אחד של 3 ביטים או יותר של 1 -> חיסור של שתי הזזות (שיטת Booth)
            // דוגמה: 14 בבינארי זה 1110 (רצף באורך 3). ניתן להציגו כ- 16 פחות 2.
            if (bitCount >= 3 && blockLength == bitCount) {
                int shiftHigher = highestSetBit + 1;
                int shiftLower = lowestSetBit;

                out.println("movq %rdi, %rax");
                out.println("shlq $" + shiftHigher + ", %rax");
                out.println("movq %rdi, %rcx");
                if (shiftLower > 0) {
                    out.println("shlq $" + shiftLower + ", %rcx");
                }
                out.println("subq %rcx, %rax");
                out.println("ret");
                return;
            }

            // פירוק לסכום של הזזות עבור כל ביט דולק, בצורה המינימלית ביותר
            boolean first = true;
            for (int i = 0; i <= highestSetBit; i++) {
                if ((k & (1 << i)) != 0) {
                    if (first) {
                        out.println("movq %rdi, %rax");
                        if (i > 0) {
                            out.println("shlq $" + i + ", %rax");
                        }
                        first = false;
                    } else {
                        out.println("movq %rdi, %rcx");
                        if (i > 0) {
                            out.println("shlq $" + i + ", %rcx");
                        }
                        out.println("addq %rcx, %rax");
                    }
                }
            }
            out.println("ret");

        } catch (IOException e) {
            System.err.println("Error writing to file: " + e.getMessage());
        }
    }
}