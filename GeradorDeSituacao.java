public class GeradorDeSituacao
{
	public static void main(String args[])
	{
		System.err.print("N de bolinhas: ");
		int n = StdIn.readInt();
		System.err.print("Semente: ");
		int seed = StdIn.readInt();
		System.err.print("Raio: ");
		double raio = StdIn.readDouble();
		System.err.print("Tempo minimo de recuperacao: ");
		double min_rec = StdIn.readDouble();
		System.err.print("Tempo maximo de recuperacao: ");
		double max_rec = StdIn.readDouble();
		System.err.print("Proporcao minima do raio de contagio: ");
		double min_spread_size = StdIn.readDouble();
		System.err.print("Proporcao maxima do raio de contagio: ");
		double max_spread_size = StdIn.readDouble();
		System.err.print("Fator de infectividade minima: ");
		double infectivity = StdIn.readDouble();
		System.err.print("Fator de infectividade geral: ");
		double base_infectivity = StdIn.readDouble();

		StdOut.println(n);
		StdOut.println(seed);
		StdOut.println(infectivity);
		StdOut.println(base_infectivity);

		StdRandom.setSeed(seed);

		for (int i = 0; i < n; i++)
		{
			StdOut.print(StdRandom.uniform(0.0, 1.0) + " ");
	        StdOut.print(StdRandom.uniform(0.0, 1.0) + " ");
	        StdOut.print(StdRandom.uniform(-0.005, 0.005) + " ");
	        StdOut.print(StdRandom.uniform(-0.005, 0.005) + " ");
	        StdOut.print(raio + " ");
	        StdOut.print(0.5 + " ");
	        if (i == 0) StdOut.print(max_rec + " ");
	        else StdOut.print(StdRandom.uniform(min_rec, max_rec) + " ");
	        StdOut.print(raio * StdRandom.uniform(min_spread_size, max_spread_size));
	        StdOut.println();
		}
		
	}
}
