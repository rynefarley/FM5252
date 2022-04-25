using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;

namespace MyApp
{
	class Program
	{
		static void Main(string[] args)
		{;

			Underlying u = new Underlying();
	
			System.Console.WriteLine("Stock Price: (ex: 97.53)");
			u.Price = Convert.ToDouble(Console.ReadLine());

			
			EuropeanOption euro = new EuropeanOption();
			System.Console.WriteLine("Enter true for call or false for put:");
			euro.IsCall =System.Convert.ToBoolean(Console.ReadLine());

				
			System.Console.WriteLine("Strike:" );
                        euro.Strike = Convert.ToDouble(Console.ReadLine());
		
			euro.Underlying = u;

			var g = euro.GetPriceAndGreeks();
			System.Console.WriteLine("Price: "+g.Price);
			System.Console.WriteLine("Delta: "+g.Delta);
			System.Console.WriteLine("Gamma: "+g.Gamma);
			System.Console.WriteLine("Vega: " +g.Vega);
			System.Console.WriteLine("Theta: " +g.Theta);
		 	System.Console.WriteLine("Rho: " +g.Rho);
			
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

		public abstract OptionResult GetPriceAndGreeks();
	}


	class OptionResult
	{
		public double Price{get;set;}
		public double Delta{get;set;}
		public double Gamma{get;set;}
		public double Vega{get;set;}
		public double Theta{get;set;}
		public double Rho{get;set;}
	}

	class Volatility
	{
		public double Vol {get;set;}
	}

	class EuropeanOption : Option
	{
		public double Strike {get;set;}
		public bool IsCall {get;set;}

		public override OptionResult GetPriceAndGreeks()
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
					
		        System.Console.WriteLine("Enter tenor: (ex 1 year =1)" );
                        p.Tenor = Convert.ToDouble(Console.ReadLine());


			System.Console.WriteLine("Enter volatility: ex(.5)" );
                        p.Volatility = Convert.ToDouble(Console.ReadLine());
			System.Console.WriteLine("Enter true for antithetic  or false for not antithetic:");
                        bool antithetic =System.Convert.ToBoolean(Console.ReadLine());



			var result = MonteCarloSimulator.GeneratePaths(p,r,antithetic);
			p.Volatility = p.Volatility+1;
			var resultVegaP  = MonteCarloSimulator.GeneratePaths(p,r,antithetic);
			p.Volatility = p.Volatility-2;
			var resultVegaN  = MonteCarloSimulator.GeneratePaths(p,r,antithetic);
			p.Volatility = p.Volatility+1;
			p.Tenor = p.Tenor*2;
			var resultTheta = MonteCarloSimulator.GeneratePaths(p,r,antithetic);
			p.Tenor=p.Tenor/2;
			p.r = p.r+1;
			var resultRho = MonteCarloSimulator.GeneratePaths(p,r,antithetic);
			p.r = p.r-1;

			
			double price = 0;
			double[] SDE = new double[Simulations];
			double deltaP = 0;
			double deltaN = 0;
			double vega = 0;
			double vegaP = 0;
			double vegaN = 0;
			double thetaP =0;
			double rhoP =0;
			
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
						
						SDE[j] = x;
						price+=x*Math.Exp(-p.r*p.Tenor)/Simulations;	
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
					if(resultVegaP.SimulatedPaths[Steps-1,j]-Strike>0)
					{
						vegaP+=(resultVegaP.SimulatedPaths[Steps-1,j]-Strike)*Math.Exp(-p.r*p.Tenor)/Simulations;
					}
				
					if(resultVegaN.SimulatedPaths[Steps-1,j]-Strike>0)
					{	vegaN+=(resultVegaN.SimulatedPaths[Steps-1,j]-Strike)*Math.Exp(-p.r*p.Tenor)/Simulations;
					}
					vega = (vegaP - vegaN) / (2*.01); 	
					if(resultTheta.SimulatedPaths[Steps-1,j]-Strike>0)
                                        {
                                                thetaP+=(resultTheta.SimulatedPaths[Steps-1,j]-Strike)*Math.Exp(-p.r*p.Tenor)/Simulations;
                                        }

					if(resultRho.SimulatedPaths[Steps-1,j]-Strike>0)
                                        {
                                                rhoP+=(resultRho.SimulatedPaths[Steps-1,j]-Strike)*Math.Exp(-p.r*p.Tenor)/Simulations;
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
					if(Strike - resultVegaP.SimulatedPaths[Steps-1,j]>0)
                                        {
                                                vegaP+=(Strike - resultVegaP.SimulatedPaths[Steps-1,j])*Math.Exp(-p.r*p.Tenor)/Simulations;
                                        }

                                        if(Strike - resultVegaN.SimulatedPaths[Steps-1,j]>0)
                                        {       vegaN+=(Strike - resultVegaN.SimulatedPaths[Steps-1,j])*Math.Exp(-p.r*p.Tenor)/Simulations;
                                        }
                                        if(Strike-resultTheta.SimulatedPaths[Steps-1,j]>0)
                                        {
                                                thetaP+=(Strike-resultVegaP.SimulatedPaths[Steps-1,j])*Math.Exp(-p.r*p.Tenor)/Simulations;
                                        }

					 if(Strike - resultRho.SimulatedPaths[Steps-1,j]>0)
                                        {
                                                rhoP+=(Strike - resultRho.SimulatedPaths[Steps-1,j])*Math.Exp(-p.r*p.Tenor)/Simulations;
                                        }


				}
			
				
			}

			deltaP = deltaP*Math.Exp(-p.r*p.Tenor);
			deltaN = deltaN*Math.Exp(-p.r*p.Tenor);
			vega = (vegaP-vegaN)/(2*.01);
			vega = vega*Math.Exp(-p.r*p.Tenor);
			double theta = thetaP-price;
			double delta = (deltaP - price)/2;
			double rho = (rhoP-price)/2;
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
                        			
			double std = 0;
			std = Math.Sqrt(error/(Simulations-1));
			double err = std/Math.Sqrt(Simulations);
			System.Console.WriteLine("Standard Error: " +err);
			return new OptionResult() {Price = price, Delta = delta, Gamma=gamma,Vega=vega,Theta=theta,Rho=rho};	
		}

	}	


	class GaussianRandoms
	{
		public double[,] NRands {get;set;}
		


		//static Tuple<double,double> BoxMuller()
        	//{	
        	//Random r = new Random();
        	//double x1 = r.NextDouble();
        	//double x2 = r.NextDouble();
        	//double z1 = Math.Sqrt(-2*Math.Log(x1)) * Math.Cos(2*Math.PI*x2);
		//double z2 = Math.Sqrt(-2*Math.Log(x1)) * Math.Sin(2*Math.PI*x2);

        	//return Tuple.Create(z1,z2);
        	//}



		public void PopulateNRands(int seed, long Steps, long Simulations )
		{
			NRands = new double[Steps,Simulations];
			Random r = new Random(seed);
		
		 for (int i = 0; i < Steps-1; i+=2)
                	for (int j = 0; j < Simulations; j ++)
                {
                    double x1, x2;
                    x1 = r.NextDouble();
                    x2 = r.NextDouble();
                    
                    NRands[i, j] = (Math.Sqrt(-2 * Math.Log(x1))) * Math.Cos(2 * Math.PI * x2);
                    NRands[i+1, j] = (Math.Sqrt(-2 * Math.Log(x1))) * Math.Sin(2 * Math.PI * x2);                 
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
		public static SimulationResult GeneratePaths(SimulationParameters p, GaussianRandoms rands, bool antithetic)
		{
			SimulationResult results = new SimulationResult();
			results.SimulatedPaths = new double[p.Steps, p.Simulations];
			if (antithetic == true)
			{
				for( int j = 0 ;j <p.Simulations-1;j+=2)
				{
				results.SimulatedPaths[0,j] = p.S0;
				results.SimulatedPaths[0,j+1] = p.S0;
					for(int i = 0;i<p.Steps-1;i++)
					{
					// results.SimulatedPaths[i+1,j] = results.SimulatedPaths[i,j] + p.r * p.Tenor/p.Steps + p.Volatility * Math.Sqrt(p.Tenor/p.Steps) * rands.NRands[i,j];
					// results.SimulatedPaths[i+1,j+1] = results.SimulatedPaths[i,j] + p.r * p.Tenor/p.Steps - p.Volatility * Math.Sqrt(p.Tenor/p.Steps) * rands.NRands[i,j];
					
					results.SimulatedPaths[i+1,j] = results.SimulatedPaths[i,j]* Math.Exp((p.r-Math.Pow(p.Volatility,2)*.5)* p.Tenor/p.Steps + p.Volatility * Math.Sqrt(p.Tenor/p.Steps) * rands.NRands[i,j]);

					results.SimulatedPaths[i+1,j+1] = results.SimulatedPaths[i,j+1]* Math.Exp((p.r-Math.Pow(p.Volatility,2)*.5)* p.Tenor/p.Steps - p.Volatility * Math.Sqrt(p.Tenor/p.Steps) * rands.NRands[i,j]);
					



                                        

					}
				}
			}
			else
			{
			for( int j = 0 ;j <p.Simulations-1;j++)
                                {
                                results.SimulatedPaths[0,j] = p.S0;
                                        for(int i = 0;i<p.Steps-1;i++)
                                        {
                                        //results.SimulatedPaths[i+1,j] = results.SimulatedPaths[i,j] + p.r * p.Tenor/p.Steps + p.Volatility * Math.Sqrt(p.Tenor/p.Steps) * rands.NRands[i,j];
                                        results.SimulatedPaths[i+1,j] = results.SimulatedPaths[i,j]* Math.Exp((p.r-Math.Pow(p.Volatility,2)*.5)* p.Tenor/p.Steps + p.Volatility * Math.Sqrt(p.Tenor/p.Steps) * rands.NRands[i,j]);

					}
                                }

			}
		
			return results;
		}	
	}

}

