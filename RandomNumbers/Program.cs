using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace RandomNumbers
{
    class Program
    {
        static void Main(string[] args)
        {
	 Console.WriteLine("Sum Twelve");
         Console.WriteLine(SumTwelve());
	 Console.WriteLine("Box-Muller");
	 Console.WriteLine(BoxMuller());
	 Console.WriteLine("Polar Rejection");
	 Console.WriteLine(Polar());
	 Console.WriteLine("Correlated Gaussian Random Numbers");
	 Console.WriteLine("Enter the correlation: (this must be between -1 and 1)");
	 double rho = Convert.ToDouble(Console.ReadLine());
	 Console.WriteLine(Correlated(rho));
	}

	
	
	static double SumTwelve()
	{
	double total = 0;
	Random r = new Random();

	for(int i = 0; i < 12;i++)
	{
	total += r.NextDouble();
	}
	
	return total - 6;
	}

	static Tuple<double,double> BoxMuller()
	{
	Random r = new Random();
	double x1 = r.NextDouble();
	double x2 = r.NextDouble();
	double z1 = Math.Sqrt(-2*Math.Log(x1)) * Math.Cos(2*Math.PI*x2);	        double z2 = Math.Sqrt(-2*Math.Log(x1)) * Math.Sin(2*Math.PI*x2);

	return Tuple.Create(z1,z2);
	}
	
	static Tuple<double,double> Polar()
	{
	Random r = new Random();
        double x1 = r.NextDouble();
        double x2 = r.NextDouble();
	double w = Math.Pow(x1,2) + Math.Pow(x2,2);
	while( w > 1)
	{
	x1 = r.NextDouble()*2-1;
        x2 = r.NextDouble()*2-1;
	w = Math.Pow(x1,2) + Math.Pow(x2,2);
	}

	double c = Math.Sqrt(-2*Math.Log(w)/w);
	double z1 = c*x1;
	double z2 = c*x2;


	return Tuple.Create(z1,z2);
	}

	static Tuple<double,double> Correlated(double rho)
	{
	Tuple<double,double> r = Polar();
	double z1 = r.Item1;
	double z2 = rho*z1 + (Math.Sqrt(1 - Math.Pow(rho,2))*r.Item2);
	if(Double.IsNaN(z2))
	{
	Console.WriteLine("That is an incorrect value. Please enter a number between -1 and 1.");

	}
	return Tuple.Create(z1,z2);
	}
    }
}
