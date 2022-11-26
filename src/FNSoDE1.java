import java.util.Scanner;

public class FNSoDE1 {

	public static void main(String[] args) {
		System.out.println("\tАнализ быстрого численного метода решения уравнения Лапласа");
		/**
		 * n - количество шагов, n1 - количество шагов для вычисления интеграла Фурье,
		 * N, T - количество уравнений в СЛАУ для уравнения теплопроводности,
		 * a - коэффициент пропорциональности, l - длина, h - высота,
		 * A, B - границы исследуемого отрезка, t - исследуемое время,
		 * Δh, Δt - шаги по пространственной и временной/температурной координатам,
		 * λ - универсальная переменная начальных условий,
		 * m - средняя квадратическая ошибка,
		 * X - массив значений пространственных координат,
		 * ΔU - вспомогательный массив значений искомой функции для записи при заданном
		 * времени, а также для вычисления средней квадратической ошибки,
		 * UT - массив значений искомой функции, вычисленных аналитическим методом,
		 * U - массив значений искомой функции, вычисленных численным методом,
		 * M - массив для составления трёхдиагональной матрицы,
		 */
		Scanner scanner = new Scanner(System.in);
		int n, Approximations = 1000, A = 0, B = 1, C = 0, D = 1;
		double ΔhX, ΔhY, m, F = 0;
		System.out.print("\nВведите размерность исследуемой области: ");
		n = scanner.nextInt() - 1;
		System.out.println("1. Решение уравнения Лапласа быстрым численным методом.\n"
				+ "2. Решение уравнения Лапласа численным методом конечных разностей.\n"
				+ "3. Выйти из программы.");
		double[] X = new double[n + 1];
		double[] Y = new double[n + 1];
		double[][] U = new double[n + 1][n + 1];
		double[][] UT = new double[n + 1][n + 1];
		double[][] UF1 = new double[n + 1][n + 1];
		double[] Yup = new double[n + 1];
		double[] Ydown = new double[n + 1];
		double[] Xleft = new double[n + 1];
		double[] Xright = new double[n + 1];
		ΔhX = (double) (B - A) / n;
		ΔhY = (double) (D - C) / n;
		X[0] = A;
		Y[0] = C;
		for (int i = 0; i < n; i++) {
			X[i + 1] = X[i] + ΔhX;
			Y[i + 1] = Y[i] + ΔhY;
		}
		for (int i = 0; i <= n; i++) {
			// Лаплас 1
			Yup[i] = Math.cos(Y[i]);
			Ydown[i] = Math.exp(1) * Math.cos(Y[i]);
			Xleft[i] = Math.exp(X[i]);
			Xright[i] = Math.exp(X[i]) * Math.cos(1);
			// Лаплас 2
			// Yup[i] = Math.sin(Math.PI * Y[i] / D);
			// Ydown[i] = 0;
			// Xleft[i] = Math.sin(Math.PI * X[i] / B);
			// Xright[i] = 0;
		}
		MatrixOfBoundaryConditions(U, Yup, Ydown, Xleft, Xright, n, A, B, C, D);
		MatrixOfBoundaryConditions(UF1, Yup, Ydown, Xleft, Xright, n, A, B, C, D);
		MatrixOfBoundaryConditions(UT, Yup, Ydown, Xleft, Xright, n, A, B, C, D);
		// Лаплас 1
		for (int i = 1; i < n; i++) {
			for (int j = 1; j < n; j++) {
				UT[i][j] = Math.exp(X[i]) * Math.cos(X[j]);
			}
		}
		// Лаплас 2
		// for (int i = 1; i < n; i++) {
		// for (int j = 1; j < n; j++) {
		// UT[i][j] = Math.sinh(Math.PI * (B - X[j]) / D) / Math.sinh(Math.PI * B / D)
		// * Math.sin(Math.PI * Y[i] / D)
		// + Math.sinh(Math.PI * (D - Y[i]) / B)
		// / Math.sinh(Math.PI * D / B) * Math.sin(Math.PI * X[j] / B);
		// }
		// }
		int K = scanner.nextInt(); // выбор значения в меню программы
		while (K != 3) { // цикл для возвращения к необходимым пунктам в меню
			if (K >= 1 && K <= 2) { // условие при вводе неверного пункта меню
				switch (K) { // объявляем список, контролируемый значением K
					case 1: // реализуется решение быстрым методом
						long millis = System.currentTimeMillis();
						m = 0;
						for (int i = 1; i < n; i++) {
							for (int j = 1; j < n; j++) {
								U[i][j] = 0;
							}
						}
						for (int i = 1; i <= n; i++) {
							for (int j = 1; j < n; j++) {
								if (i != n)
									UF1[i][j] = (UF1[i - 1][j] + UF1[i - 1][j - 1] + UF1[i][j - 1]) / 3;
								if (i != 1)
									U[i - 1][j] += (UF1[i - 2][j] + UF1[i][j] + UF1[i - 1][j - 1] + UF1[i - 1][j + 1]) / 16;
							}
						}
						for (int i = 1; i <= n; i++) {
							for (int j = n - 1; j >= 1; j--) {
								if (i != n)
									UF1[i][j] = (UF1[i - 1][j] + UF1[i - 1][j + 1] + UF1[i][j + 1]) / 3;
								if (i != 1)
									U[i - 1][j] += (UF1[i - 2][j] + UF1[i][j] + UF1[i - 1][j - 1] + UF1[i - 1][j + 1]) / 16;
							}
						}
						for (int i = n - 1; i >= 0; i--) {
							for (int j = 1; j < n; j++) {
								if (i != 0)
									UF1[i][j] = (UF1[i + 1][j] + UF1[i + 1][j - 1] + UF1[i][j - 1]) / 3;
								if (i != n - 1)
									U[i + 1][j] += (UF1[i][j] + UF1[i + 2][j] + UF1[i + 1][j - 1] + UF1[i + 1][j + 1]) / 16;
							}
						}
						for (int i = n - 1; i >= 0; i--) {
							for (int j = n - 1; j >= 1; j--) {
								if (i != 0)
									UF1[i][j] = (UF1[i + 1][j] + UF1[i + 1][j + 1] + UF1[i][j + 1]) / 3;
								if (i != n - 1)
									U[i + 1][j] += (UF1[i][j] + UF1[i + 2][j] + UF1[i + 1][j - 1] + UF1[i + 1][j + 1]) / 16;
							}
						}
						for (int k = 1; k <= Approximations; k++) {
							for (int i = 1; i < n; i++) {
								for (int j = 1; j < n; j++) {
									if (k % 2 != 0) {
										UF1[i][j] = (U[i - 1][j] + U[i + 1][j] + U[i][j - 1] + U[i][j + 1]) / 4;
									} else {
										U[i][j] = (UF1[i - 1][j] + UF1[i + 1][j] + UF1[i][j - 1] + UF1[i][j + 1]) / 4;
									}
								}
							}
						}
						if (Approximations % 2 != 0) {
							U = UF1;
						}
						System.out.println("Решение уравнения Лапласа быстрым численным методом:");
						// PrintTable(U, X, Y, n);
						System.out.println("Точное решение уравнения Лапласа:");
						// PrintTable(UT, X, Y, n);
						for (int i = 1; i < n; i++) {
							for (int j = 1; j < n; j++) {
								m += Math.abs((UT[i][j] - U[i][j]) / UT[i][j] * 100);
							}
						}
						m = m / Math.pow((n - 1), 2);
						System.out.println("Погрешность вычислений в процентах: " + String.format("%.6f", m));
						System.out
								.println("Время проведённых вычислений: " + ((System.currentTimeMillis() - millis)) + " мс.");
						break;
					case 2:
						// millis = System.currentTimeMillis();
						m = 0;
						FiniteDifferenceMethodLP(U, Xleft, Xright, Yup, Ydown, ΔhX, ΔhY, F, n);
						System.out.println("Решение уравнения Лапласа численным методом конечных разностей:");
						// PrintTable(U, X, Y, n);
						System.out.println("Точное решение уравнения Лапласа:");
						// PrintTable(UT, X, Y, n);
						for (int i = 1; i < n; i++) {
							for (int j = 1; j < n; j++) {
								m += Math.abs((UT[i][j] - U[i][j]) / UT[i][j] * 100);
							}
						}
						m = m / Math.pow((n - 1), 2);
						System.out.println("Погрешность вычислений в процентах: " + String.format("%.6f", m));
						// System.out
						// .println("Время проведённых вычислений: " + ((System.currentTimeMillis() -
						// millis)) + " мс.");
						break;
				}
			} else {
				System.out.println("Введите значение от 1 до 3!");
			}
			K = scanner.nextInt();
		}
		System.out.println("Программа завершена.");
		scanner.close();
	}

	public static void FiniteDifferenceMethodLP(double U[][], double[] Xleft, double[] Xright,
			double[] Yup, double[] Ydown, double F, double hx, double hy, int n) { // метод конечных разностей
		int p = (n - 1) * (n - 1);
		double[][] M = new double[p][p + 1];
		double[] ΔU = new double[p];
		for (int i = 0; i < n - 1; i++) {
			for (int j = i * (n - 1); j < (i + 1) * (n - 1); j++) {
				M[j][j] = 4;
				if (j % (n - 1) != n - 2) {
					M[j + 1][j] = -1;
					M[j][j + 1] = -1;
				}
				if (i != (n - 2)) {
					M[j][j + n - 1] = -1;
				}
				if (i != 0) {
					M[j][j - n + 1] = -1;
				}
			}
		}
		M[0][p] = Xleft[1] + Yup[1];
		M[n - 2][p] = Xright[1] + Yup[n - 1];
		M[p - n + 1][p] = Xleft[n - 1] + Ydown[1];
		M[p - 1][p] = Xright[n - 1] + Ydown[n - 1];
		for (int i = 1; i < n - 2; i++) {
			M[i][p] = Yup[i + 1];
			M[i * (n - 1)][p] = Xleft[i + 1];
			M[i * (n - 1) + (n - 2)][p] = Xright[i + 1];
			M[i + p - n + 1][p] = Ydown[i + 1];
		}
		if (F != 0) {
			for (int i = 0; i < p; i++) {
				M[i][p] -= hx * hy * F;
			}
		}
		methodGauss(M, p, ΔU);
		for (int i = 1; i < n; i++) {
			for (int j = 1; j < n; j++) {
				U[i][j] = ΔU[j - 1 + (i - 1) * (n - 1)];
			}
		}
	}

	public static void MatrixOfBoundaryConditions(double U[][], double Yup[], double Ydown[], double Xleft[],
			double Xright[], int n, int A, int B, int C, int D) { // присвоение таблицам граничных значений
		for (int i = 0; i <= n; i++) {
			U[C][i] = Yup[i];
			U[D * n][i] = Ydown[i];
			U[i][A] = Xleft[i];
			U[i][B * n] = Xright[i];
		}
	}

	public static void PrintTable(double U[][], double X[], double Y[], int n) { // функция для вывода таблицы значений
		System.out.print("╔═══════╦");
		for (int i = 0; i < n; i++) {
			System.out.print("══════════╦");
		}
		System.out.print("══════════╗\n║ " + "U(x,y)" + "");
		for (int i = 0; i < n + 1; i++) { // цикл, заполняющий шапку таблицы
			System.out.print("║   " + String.format("%.2f", X[i]) + "   ");
		}
		System.out.print("║\n╠═══════╬");
		for (int i = 0; i < n; i++) {
			System.out.print("══════════╬");
		}
		System.out.print("══════════╣\n");
		for (int i = 0; i < n + 1; i++) { // цикл, проходящий строки таблицы
			System.out.print("║ " + String.format("%.2f", Y[i]) + "\t");
			for (int j = 0; j < n + 1; j++) { // цикл, заполняющий ячейки текущей строки
				System.out.print("║ " + String.format("%.6f", U[i][j]) + " ");
			}
			System.out.println("║");
		}
		System.out.print("╚═══════╩");
		for (int i = 0; i < n; i++) {
			System.out.print("══════════╩");
		}
		System.out.print("══════════╝\n");
	}

	public static void methodGauss(double A[][], int n, double X[]) { // функция, реализующая метод Гаусса
		double[] B = new double[n];
		for (int i = 0; i < n; i++) {
			B[i] = A[i][n];
		}
		// прямой ход (или приведение матрицы к треугольному виду)
		for (int p = 0; p < n; p++) { // с помощью цикла перебираем все уравнения
			int max = p; // задаём максимальное значение первого ненулевого коэффициента в уравнении
			for (int i = p + 1; i < n; i++) { // пробегаем по столбцам и ищем максимальное значение
				if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
					max = i;
				}
			}
			double[] temp = A[p];
			A[p] = A[max]; // первой строке СЛАУ присваиваем строку с максимальным значением первого
								// ненулевого коэффициента
			A[max] = temp;
			double t = B[p]; // ставим на место выпавшее уравнение
			B[p] = B[max];
			B[max] = t;
			for (int i = p + 1; i < n; i++) { // цикл, преобразующий строку так, что её первый элемент, считая над
															// первым
															// ненулевым элементом, становится равен нулю
				double alpha = A[i][p] / A[p][p]; // получение коэффициента, на который будет умножаться строка,
																// отнимаемая
																// от всех последующих
				B[i] -= alpha * B[p]; // отдельно производим вычитание для свободных членов
				for (int j = p; j < n; j++) {
					A[i][j] -= alpha * A[p][j]; // производим вычитание коэффициента, умноженного на строку так, что её
															// первый элемент, считая над первым ненулевым элементом, становится
															// равен
															// нулю
				}
			}
		}
		// обратный ход
		for (int i = n - 1; i >= 0; i--) { // обратный ход осуществим через два цикла; пройдём все строки СЛАУ в обратном
														// порядке
			double sum = 0.0; // вводим счётчик
			for (int j = i + 1; j < n; j++) { // цикл заключается в получении в счётчик значения, при котором при
															// подстановке на одну строчку вверх нашего найденного значения b[i] мы
															// будем получать значение, на которое нужно отнять свободный член и
															// поделить это на свободный коэффициент в строке, чтобы получить следующее
															// значение b[i-1]
				sum += A[i][j] * X[j];
			}
			X[i] = (B[i] - sum) / A[i][i]; // получаем значение b[i], то есть делим свободный член (отнимать счётчик при
														// последующих манипуляциях) на коэффициент при b[i]
		}
	}
}