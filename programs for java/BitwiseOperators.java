public class BitwiseOperators {
    public static void main(String[] args) {
        int a = 6; //binary: 0110
        int b = 4; //binary: 0100
        System.out.println("a=" + a);
        System.out.println("b=" + b);
        
        System.out.println("a&b=" + (a & b)); // bitwise AND
        System.out.println("a|b=" + (a | b)); //bitwise OR
        System.out.println("a^b=" + (a ^ b)); // bitwise XOR
        System.out.println("~a=" +(~a)); // bitwise NOT
        System.out.println("a<<1=" + (a << 1)); //bitwise Left shift
        System.out.println("a>>1=" + (a >>1));


    }
}
