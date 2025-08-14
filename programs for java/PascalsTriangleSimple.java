public class PascalsTriangleSimple {

    public static void main(String[] args) {
        int n = 5; // Number of rows
        int[][] triangle = new int[n][n];

        // Build triangle using simple logic
        for (int i = 0; i < n; i++) {
            // Each row starts and ends with 1
            triangle[i][0] = 1;
            triangle[i][i] = 1;

            // Fill inside values: sum of two elements above
            for (int j = 1; j < i; j++) {
                triangle[i][j] = triangle[i - 1][j - 1] + triangle[i - 1][j];
            }
        }

        // Print the triangle
        for (int i = 0; i < n; i++) {
            // Print spaces for formatting
            for (int k = 0; k < n - i; k++) {
                System.out.print(" ");
            }

            // Print row elements
            for (int j = 0; j <= i; j++) {
                System.out.print(triangle[i][j] + " ");
            }
            System.out.println();
        }
    }
}
