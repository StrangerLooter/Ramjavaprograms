import java.util.Scanner;

public class SumNatural {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);

        // Ask user for input before reading it
        System.out.println("Enter the positive number:");
        int n = sc.nextInt();

        if (n <= 0) {
            System.out.println("Please enter a positive number");
        } else {
            // Using formula
            int sum = n * (n + 1) / 2;
            System.out.println("The sum of first " + n + " natural numbers is: " + sum);
        }
    }
}
