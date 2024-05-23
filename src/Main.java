
import java.math.*;//TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
import java.text.DecimalFormat;

// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
public class Main {

    public static  Vector solveR(Matrix R, Vector b, int flag)
    {
        double sum = 0;
        Vector x = new Vector(b.getSize());
    if(flag == 0)
        for(int i = 0;i<b.getSize();++i)
        {
            for(int j = 0;j<i;++j)
            {
                sum+= R.getElemnt(i,j)*x.getElement(j);
            }
            sum = b.getElement(i) - sum;
            sum = sum / R.getElemnt(i,i);
            x.setElement(i,sum);
            sum = 0;
        }
    if(flag == 1)
        for(int i = b.getSize()-1;i>=0;--i)
        {
            for(int j = b.getSize()-1;j>i;--j)
            {
                sum+= R.getElemnt(i,j)*x.getElement(j);
            }
            sum = b.getElement(i) - sum;
            sum = sum / R.getElemnt(i,i);
            x.setElement(i,sum);
            sum = 0;
        }
        return x;
    }
    public static Vector yacobi(Matrix A, double error)
    {
        Matrix D = new Matrix(A.getSize());
        Matrix Q = new Matrix(A.getSize(),1);
        Vector eigenvalyes = new Vector(A.getSize());
        double tan, sin, cos,tau;
        D.equal(A);
        while (true){
        for(int i = 0;i<D.getSize();++i)
        {
            for(int j = 0;j<i;++j)
            {
                if((D.getElemnt(i,i) - D.getElemnt(i,i)*D.delta() - D.delta()) < error) return D.trace();
                if(D.getElemnt(i,j)!= 0.0)
                {
                    tau = (D.getElemnt(i,i) - D.getElemnt(j,j))/D.getElemnt(i,j);
                    tan = 1/(tau + Math.signum(tau)*Math.sqrt(tau*tau + 1));
                    cos = 1/(Math.sqrt(tan*tan +1));
                    sin = tan * cos;
                    Q.setElement(i,i,cos);
                    Q.setElement(j,j,cos);
                    Q.setElement(i,j,sin);
                    Q.setElement(j,i,-sin);
                    D.equal(Q.transposition().multMatrix(D.multMatrix(Q)));

                }

            }
        }
    }
    }
    public static Vector QR(Matrix A, double error)
    {
        Matrix B1 = new Matrix(A.getSize(),0);
        Matrix B = new Matrix(A.getSize());
        Vector b = new Vector(A.getSize());
        B1.equal(A);
        while (true)
        {
            B.equal(rollingMethod11(B1,b));
            if(B1.trace().sumVector(B.trace().multOnScal(-1)).chebNorm()<error)
            {
                return B.trace();
            }
            B1.equal(B);
        }
    }
    public static Vector conjugateGradient(Matrix A,Vector b,double epsilon)
    {
        int k = 0;
        double lamda = Math.pow(yacobi(A,epsilon).chebNorm(),-1),step,v;
        Vector x = new Vector(A.getSize(),0);
        Vector discrepancy = b.sumVector(A.multOnVector(x).multOnScal(-1.0));
        Vector discrepancy1 = new Vector(b.getSize());
        discrepancy1.equal(discrepancy);
        Vector s = new Vector(b.getSize());
        Vector g = new Vector(b.getSize());
        s.equal(discrepancy);
        while (lamda*(A.multOnVector(x).sumVector(b.multOnScal(-1))).norm()>epsilon)
        {
            g.equal(A.multOnVector(s));
            step = discrepancy.scalMult(discrepancy)/s.scalMult(g);
            x.equal(x.sumVector(s.multOnScal(step)));
            discrepancy.equal(discrepancy.sumVector(g.multOnScal(step*(-1))));
            v = discrepancy.scalMult(discrepancy)/discrepancy1.scalMult(discrepancy1);
            s.equal(discrepancy.sumVector(s.multOnScal(v)));
            discrepancy1.equal(discrepancy);
            k++;

        }
        System.out.println(k);
        return x;

    }
    public static Vector fallDawn(Matrix A,Vector b,double epsilon)
    {
        double step,lamda = Math.pow(yacobi(A,epsilon).chebNorm(),-1);
        Vector x = new Vector(b.getSize(),1);
        Vector discrepancy;
        int numberStep =0;

        //(0.5*A.multOnVector(x).scalMult(x) - b.scalMult(x))<=epsilon

        while (lamda*(A.multOnVector(x).sumVector(b.multOnScal(-1))).norm()>epsilon)
        {
            discrepancy = A.multOnVector(x).sumVector(b.multOnScal(-1.0));
            step = (Math.pow(discrepancy.norm(),2)/A.multOnVector(discrepancy).scalMult(discrepancy));
            x.equal(x.sumVector(discrepancy.multOnScal(-step)));
            numberStep++;
        }
        return x;
    }
    public static Matrix rollingMethod11(Matrix A,Vector b)
    {
        Matrix Q = new Matrix(A.getSize(),1);
        Matrix R = new Matrix(A.getSize(),1);
        Matrix roll;
        R.equal(A);
        double cosa,sina;
        for(int i = 0;i<A.getSize();++i)
        {
            for (int j = i+1;j<A.getSize();++j)
            {
                roll = new Matrix(A.getSize(), 1);
                cosa = R.getElemnt(i,i)/(Math.sqrt(R.getElemnt(i,i)*R.getElemnt(i,i) +
                        R.getElemnt(j,i)*R.getElemnt(j,i)));
                sina = (-R.getElemnt(j,i))/(Math.sqrt((R.getElemnt(i,i)*R.getElemnt(i,i) +
                        R.getElemnt(j,i)*R.getElemnt(j,i))));
                roll.setElement(i,i,cosa);
                roll.setElement(j,i,sina);
                roll.setElement(i,j,-sina);
                roll.setElement(j,j,cosa);
                R.equal(roll.multMatrix(R));
                Q.equal(roll.multMatrix(Q));

            }
        }
        return R.multMatrix(Q.transposition());
    }
    public static Vector rollingMethod(Matrix A,Vector b, Matrix B)
    {
        Matrix Q = new Matrix(A.getSize(),1);
        Matrix R = new Matrix(A.getSize(),1);
        Matrix roll;
        R.equal(A);
        double cosa,sina;
        for(int i = 0;i<A.getSize();++i)
        {
         for (int j = i+1;j<A.getSize();++j)
         {
             roll = new Matrix(A.getSize(), 1);
             cosa = R.getElemnt(i,i)/(Math.sqrt(R.getElemnt(i,i)*R.getElemnt(i,i) +
                                            R.getElemnt(j,i)*R.getElemnt(j,i)));
             sina = (-R.getElemnt(j,i))/(Math.sqrt((R.getElemnt(i,i)*R.getElemnt(i,i) +
                                            R.getElemnt(j,i)*R.getElemnt(j,i))));
             roll.setElement(i,i,cosa);
             roll.setElement(j,i,sina);
             roll.setElement(i,j,-sina);
             roll.setElement(j,j,cosa);
             R.equal(roll.multMatrix(R));
             Q.equal(roll.multMatrix(Q));

         }
        }
        System.out.println(R);
        System.out.println(Q);
        B = R.multMatrix(Q);
        B.equal(R.multMatrix(Q));
        System.out.println("//////////////////////");
        return solveR(R,Q.multOnVector(b),1);
    }
    public static double rollingMethodmin(Matrix A)
    {
        Matrix R = new Matrix(A.getSize(),1);
        Matrix roll;
        R.equal(A);
        double cosa,sina,min = 0;
        for(int i = 0;i<A.getSize();++i)
        {
            for (int j = i+1;j<A.getSize();++j)
            {
                roll = new Matrix(A.getSize(), 1);
                cosa = R.getElemnt(i,i)/(Math.sqrt(R.getElemnt(i,i)*R.getElemnt(i,i) +
                        R.getElemnt(j,i)*R.getElemnt(j,i)));
                sina = (-R.getElemnt(j,i))/(Math.sqrt((R.getElemnt(i,i)*R.getElemnt(i,i) +
                        R.getElemnt(j,i)*R.getElemnt(j,i))));
                roll.setElement(i,i,cosa);
                roll.setElement(j,i,sina);
                roll.setElement(i,j,-sina);
                roll.setElement(j,j,cosa);
                R.equal(roll.multMatrix(R));

            }
        }
        min = R.getElemnt(0,0);
        for(int i = 0;i< R.getSize();++i)
        {
            if(min< R.getElemnt(i,i))
            min = R.getElemnt(i,i);
        }
        return Math.abs(min);
    }
    public static Vector holess(Matrix A, Vector b)
    {
        double sum = 0;
        Vector y;
        Vector x;
        Matrix hol = new Matrix(A.getSize(),0);

        for(int j = 0;j<hol.getSize();j++)
        {
            sum = 0;
            for(int k = 0;k<j;++k)
            {
                sum += hol.getElemnt(j,k) * hol.getElemnt(j,k);
            }
            sum = A.getElemnt(j,j) - sum;
            if(sum<0)
            {
                System.out.printf("Попытка извлечь корень из отрицательного %.16e\n", sum);
                break;
            }
            sum = Math.sqrt((sum));
            hol.setElement(j,j,sum);
            sum = 0;
            for(int i = j+1;i< A.getSize();i++)
            {
                for(int k = 0;k<j;++k)
                {
                    sum += hol.getElemnt(i,k) * hol.getElemnt(j,k);
                }
                sum = A.getElemnt(i,j) - sum;
                sum = sum / hol.getElemnt(j,j);
                hol.setElement(i,j,sum);
                sum = 0;

            }
        }
        y = solveR(hol,b,0);
        x = solveR(hol.transposition(),y,1);
        return x;
    }
    public static void main(String[] args) {
        int n =5, k = 13;
        DecimalFormat df = new DecimalFormat("0.000000000000000E0");
        double iterations = 0.0000001;
        double time,timeall=0;
        double[] vector_massive = new double[n];
        for(int i = 0;i<n;++i)
        {
            vector_massive[i] = Math.pow(-1.0,i);
        }
        Matrix Gilbert = new Matrix(n,2);
        Matrix Gilbert1 = new Matrix(n,3);
        Matrix A = new Matrix(n,2);
        Vector answers = new Vector(vector_massive,n);
//
//        System.out.println(Gilbert);
//        System.out.println(answers);
//        //System.out.println(Gilbert.determinant(Gilbert1,n));
//        time = System.currentTimeMillis();
//        System.out.println("Holess\n");
//        System.out.println("Answer = " + holess(Gilbert1,answers));
//        System.out.println("Left part = " + Gilbert.multOnVector(holess(Gilbert1,answers)));
//        System.out.println("Run time = " + (System.currentTimeMillis() - time)/2 + "\n");
//        System.out.println("Nevyazka = " + Gilbert.multOnVector(holess(Gilbert1,answers)).sumVector(answers.multOnScal(-1.0)));
//        System.out.println("Norm of nevyazka = " + (Gilbert.multOnVector(holess(Gilbert,answers)).sumVector(answers.multOnScal(-1.0))).norm());
//        System.out.println("\n\n _____________________________ \n\n");
//        System.out.println("Rolling method\n");
//        time = System.currentTimeMillis();
//        System.out.println("Answer = " + rollingMethod(Gilbert1,answers));
//        System.out.println("Left part = " + Gilbert.multOnVector(rollingMethod(Gilbert1,answers)));
//        System.out.println("Run time = " + (System.currentTimeMillis() - time)/2);
//        System.out.println("Nevyazka = " + Gilbert.multOnVector(rollingMethod(Gilbert1,answers)).sumVector(answers.multOnScal(-1.0)));
//        System.out.printf("Norm of nevyazka = %.16e\n", (Gilbert.multOnVector(rollingMethod(Gilbert1,answers)).sumVector(answers.multOnScal(-1.0))).norm());
//        System.out.printf("Разница = %.16e\n", (Gilbert.multOnVector(rollingMethod(Gilbert1,answers)).sumVector(answers.multOnScal(-1.0)).norm()
//                - Gilbert.multOnVector(holess(Gilbert,answers)).sumVector(answers.multOnScal(-1.0)).norm()));
////        System.out.println("Falling Dawn method\n");
//        System.out.println("Answer = " + fallDawn(Gilbert1,answers,iterations));
//        System.out.println("Left part = " + Gilbert.multOnVector(fallDawn(Gilbert1,answers,iterations)));
//        System.out.println("Run time = " + (System.currentTimeMillis() - time));
//        System.out.println("Nevyazka = " + Gilbert.multOnVector(fallDawn(Gilbert1,answers,iterations)).sumVector(answers.multOnScal(-1.0)));
//        System.out.println("Norm of nevyazka = " + (Gilbert.multOnVector(fallDawn(Gilbert1,answers,iterations)).sumVector(answers.multOnScal(-1.0))).norm());
//        System.out.println("\n Разница = " + (Gilbert.multOnVector(fallDawn(Gilbert1,answers,iterations)).sumVector(answers.multOnScal(-1.0)).norm()
//                - Gilbert.multOnVector(holess(Gilbert,answers)).sumVector(answers.multOnScal(-1.0)).norm()));
        time = 0;
//        System.currentTimeMillis();
//        for(int i = 0;i<k;++i)
//        {
//            time = System.currentTimeMillis();
//             Vector c = holess(Gilbert, answers);
//             time = System.currentTimeMillis() - time;
//             timeall += time;
//        }
//        System.out.printf("Run time Holess = %.16e\n", timeall/k);
//        timeall =0;
//        for(int i = 0;i<k;++i)
//        {
//            time = System.currentTimeMillis();
//            Vector c = rollingMethod(Gilbert, answers);
//            time = System.currentTimeMillis() - time;
//            timeall += time;
//        }
//        System.out.printf("Run time Roll = %.16e\n", timeall/k);

        System.out.println(A);
        for(int i = 0;i<100;++i)
        {
            A.equal(rollingMethod11(A,answers));
        }
        System.out.println(A);
        System.out.println(yacobi(A,iterations));
        System.out.println(QR(A,iterations));
        System.out.println(yacobi(A,iterations).sumVector(QR(A,iterations).multOnScal(-1)).norm());
        System.out.println("***************************************************************");
        System.out.println(holess(A,answers));
        System.out.println(fallDawn(A,answers,iterations));
        System.out.println(conjugateGradient(A,answers,iterations));



    }
}