import java.util.Scanner;
public class VotingEligibility {
    public static void main(String[] args){
        Scanner sc=new Scanner(System.in);
        
        System.out.println("enter the age of person: ");
         int age = sc.nextInt();
         if(age>=18){
            System.out.println("you are eligible for voting");
         }
         else {
             System.out.println("not eligible");
         }
         sc.close();
         }

    }
    

