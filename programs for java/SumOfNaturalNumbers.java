import java.util.Scanner;

public class SumOfNaturalNumbers {
    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);

        // Ask user for input
        System.out.print("Enter a positive integer: ");
        int n = input.nextInt();

        // Check if input is valid
        if (n <= 0) {
            System.out.println("Please enter a positive integer.");
        } else {
            // Using formula: Sum = n * (n + 1) / 2
            int sum = n * (n + 1) / 2;

            // Display the result
            System.out.println("The sum of first " + n + " natural numbers is: " + sum);
        }

        input.close();
    }
}
