using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;

namespace MyApp
{
	class Program
	{
		static void Main(string[] args)
		{
			Exchange e = new Exchange();
			e.Name = "New York Stock Exchange";
			e.Symbol = "NYSE";

			Underlying u = new Underlying();
			u.Exchange = e;
			u.Symbol = "GE";
			u.Price = 100;
			
			EuropeanOption euro = new EuropeanOption();
			euro.IsCall=true;
			euro.Strike=100;
			euro.ExpirationDate = new DateTime(2022,06,10);
			euro.Underlying = u;

			Volatility vol = new Volatility() {Vol = 50};

			var g = euro.GetPriceAndGreeks(vol);
			System.Console.WriteLine(g.Price);
			System.Console.WriteLine(g.Delta);
			System.Console.WriteLine(g.Gamma);
						
		}
	}

	class RatePoint
	{
		public double Tenor {get; set;}
		public double Rate {get; set;}
	}

	class YieldCurve
	{
		public List<RatePoint> CurvePoints {get;set;}
	}
	
	class Exchange
	{
		public string Name { get; set;}
		public string Symbol { get; set;}
	}

	class Underlying 
	{
		public string Symbol {get; set;}
		public Exchange Exchange {get; set;}
	        public double Price {get;set;}	
	}

	abstract class Option
	{
		public DateTime ExpirationDate {get;set;}
		public Underlying Underlying {get;set;}

		public abstract OptionResult GetPriceAndGreeks(Volatility vol);
	}


	class OptionResult
	{
		public double Price{get;set;}
		public double Delta{get;set;}
		public double Gamma{get;set;}
	}

	class Volatility
	{
		public double Vol {get;set;}
	}

	class EuropeanOption : Option
	{
		public double Strike {get;set;}
		public bool IsCall {get;set;}

		public override OptionResult GetPriceAndGreeks( Volatility vol)
		{
			System.Console.WriteLine("Enter number of steps");
			long Steps = long.Parse(Console.ReadLine());
			System.Console.WriteLine("Enter number of simulations:");
			long Simulations = long.Parse(Console.ReadLine());
			GaussianRandoms r = new GaussianRandoms();
			r.PopulateNRands(1,Steps,Simulations);

			SimulationParameters p = new SimulationParameters();
			p.S0 = Underlying.Price;
			System.Console.WriteLine("Enter the risk free rate: (ex: 0.05)");
			p.r = Convert.ToDouble(Console.ReadLine());
			p.Steps = Steps;
			p.Simulations = Simulations; 
			p.Tenor = (ExpirationDate - DateTime.Today).Days / 365d;
			p.Volatility = vol.Vol; 

			var result = MonteCarloSimulator.GeneratePaths(p,r);
			
			double price = 0;
			double[] SDE = new double[Simulations];
			double deltaP = 0;
			double deltaN = 0;
			for(int j = 0; j<Simulations; j++)
			{
				if(IsCall== true)
				{
					if(j==0)
					{
						System.Console.WriteLine("European Call Value:");
					}
					if(result.SimulatedPaths[Steps-1,j]-Strike > 0)
					{	
						double x=result.SimulatedPaths[Steps-1,j]-Strike;
						price += x/Simulations;
						SDE[j] = x;
						
					}
					else
					{
						SDE[j] = 0;
					}
                                        if(result.SimulatedPaths[Steps-1,j]+2-Strike >0)
                                        {
                                                deltaP +=(result.SimulatedPaths[Steps-1,j]+2-Strike)/Simulations;
                                        }
                                        if(result.SimulatedPaths[Steps-1,j]-2-Strike>0)
                                        {
                                                deltaN += (result.SimulatedPaths[Steps-1,j]-2-Strike)/Simulations;
                                        }

				}
				else if(IsCall== false)
				{
					if(j==0)
					{
						System.Console.WriteLine("European Put Value:");
					}
					if(Strike - result.SimulatedPaths[Steps-1,j] > 0)
					{
						double x = Strike - result.SimulatedPaths[Steps-1,j];
				       		price += x/Simulations;
						SDE[j] = x;	
					}
					else
					{
						SDE[j] = 0;
					}
					if(Strike - result.SimulatedPaths[Steps-1,j] +2 > 0)
					{
						deltaP +=(Strike - result.SimulatedPaths[Steps-1,j]+2)/Simulations;
					}
					if(Strike - result.SimulatedPaths[Steps-1,j]-2 > 0)
					{
						deltaN += (Strike - result.SimulatedPaths[Steps-1,j]-2)/Simulations;
					}
				}
			
			}

			deltaP = deltaP*Math.Exp(-p.r*p.Tenor);
			deltaN = deltaN*Math.Exp(-p.r*p.Tenor);
			price = price*Math.Exp(-p.r*p.Tenor);
			double delta = (deltaP - price)/2;
			if(IsCall == false)
			{
				delta = (deltaN - price)/2;
	
			}
			double gamma = (deltaP - 2*price + deltaN)/4;
			double error = 0;
			for(int j=0;j<Simulations;j++)
			{
				error += Math.Pow(SDE[j]*Math.Exp(-p.r*p.Tenor) - price,2);
			}
                        			
			double err = Math.Sqrt(error/(Simulations-1));
			System.Console.WriteLine("Standard Error: " +err);
			return new OptionResult() {Price = price, Delta = delta, Gamma=gamma};	
		}

	}	


	class GaussianRandoms
	{
		public double[,] NRands {get;set;}
		


		static Tuple<double,double> BoxMuller()
        	{	
        	Random r = new Random();
        	double x1 = r.NextDouble();
        	double x2 = r.NextDouble();
        	double z1 = Math.Sqrt(-2*Math.Log(x1)) * Math.Cos(2*Math.PI*x2);
		double z2 = Math.Sqrt(-2*Math.Log(x1)) * Math.Sin(2*Math.PI*x2);

        	return Tuple.Create(z1,z2);
        	}



		public void PopulateNRands(int seed, long rows, long cols )
		{
			NRands = new double[rows+2,cols];
			Random r = new Random(seed);

			for(int i = 0; i < rows/2+1; i++)
			{
				for(int j = 0; j < cols; j++)
				{	
					Tuple<double,double> a = BoxMuller();
			                NRands[2*i,j] = a.Item1;
                			NRands[2*i+1,j] = a.Item2;
				}
			}
		}
	
	}

	class SimulationResult
	{
		public double[,] SimulatedPaths {get;set;}
	}
	
	class SimulationParameters
	{
		public double S0 {get;set;}
		public double r {get;set;}
		public long Steps {get;set;}
		public long Simulations {get;set;}
		public double Tenor {get;set;}
		public double Volatility {get;set;}
	}

	static class MonteCarloSimulator
	{
		public static SimulationResult GeneratePaths(SimulationParameters p, GaussianRandoms rands)
		{
			SimulationResult results = new SimulationResult();
			results.SimulatedPaths = new double[p.Steps, p.Simulations];
			for( int j = 0 ;j <p.Simulations-1;j++)
			{
			results.SimulatedPaths[0,j] = p.S0;
				for(int i = 0;i<p.Steps-1;i++)
				{
					results.SimulatedPaths[i+1,j] = results.SimulatedPaths[i,j] + p.r * p.Tenor/p.Steps + p.Volatility * Math.Sqrt(p.Tenor/p.Steps) * rands.NRands[i,j];
				}
			}
		
			return results;
		}	
	}

}

