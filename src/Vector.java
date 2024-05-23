import java.util.Arrays;
import java.math.*;
import java.text.DecimalFormat;
import java.util.Locale;
public class Vector {
    private double[] vector;
    private int size;
    public Vector(int userSize)
    {
        vector = new double[userSize];
        size = userSize;
    }
    public Vector(double[] a,int size)
    {
        vector = new double[size];
        for(int i = 0;i<size;++i)
        {
            vector[i] = a[i];
        }
        this.size = size;
    }
    public Vector(int size, double element)
    {
        vector = new double[size];
        for(int i = 0;i<size;++i)
        {
            vector[i] = element;
        }
        this.size = size;
    }

    public double getElement(int number) {
        return vector[number];
    }

    public void setElement(int number,double element) {
        this.vector[number] = element;
    }
    public double norm()
    {
        double sum = 0;
        for(int i = 0;i<size;++i)
        {
            sum+=vector[i]*vector[i];
        }
        return Math.sqrt(sum);
    }

    public int getSize() {
        return size;
    }
    public Vector multOnScal(double scal)
    {
        Vector v = new Vector(this.getSize());
        for(int i = 0;i<size;++i)
        {
            v.setElement(i,this.getElement(i)*scal);
        }
        return v;
    }
    public Vector multOnVector(Vector b)
    {
        Vector answer = new Vector(this.getSize());
        for(int i = 0;i<this.getSize();++i)
        {
            answer.setElement(i,this.getElement(i)* b.getElement(i));
        }
        return answer;
    }
    public Vector sumVector(Vector u)
    {
      Vector answer = new Vector(this.getSize());
      for (int i =0;i<this.getSize();++i)
      {
          answer.setElement(i, this.getElement(i) + u.getElement(i));
      }
      return answer;
    }
    public double scalMult(Vector u)
    {
        return this.multOnVector(u).norm();
    }
    public void equal(Vector b)
    {
        for(int i = 0; i< this.getSize();++i)
        {
            this.setElement(i,b.getElement(i));
        }
    }
    public double chebNorm()
    {
        double max = 0;
        for (int i = 0;i<this.getSize();++i)
        {
            if(max<Math.abs(this.getElement(i)))
            {
                max = Math.abs(this.getElement(i));
            }
        }
        return max;
    }


    @Override
    public String toString() {
        String str = new String();
        DecimalFormat df = new DecimalFormat("0.000000000000000E0");
        str ="Vector{";
        for(int i = 0;i<size;++i)
            str = str  + df.format(vector[i])+"  ";

        str += "}\n";
        return str;
    }
}

