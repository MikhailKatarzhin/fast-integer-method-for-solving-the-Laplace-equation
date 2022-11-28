import java.util.Scanner;

public class FNSoDE1 {

	static Scanner scanner = new Scanner(System.in);
	static double calculationError;
	static int approximations;
	static final int DEFAULT_TABLE_SIZE;
	static int tableSize;
	static double[][] U;
	static double[][] UT;
	static double[][] UF1;
	static double[] Yup;
	static double[] Ydown;
	static double[] Xleft;
	static double[] Xright;
	static double ΔhX;
	static double ΔhY;
	static double F;
	static int A = 0, B = 1, C = 0, D = 1;
	static boolean printingTable = false;

	static {
		approximations = 1000;
		DEFAULT_TABLE_SIZE = 3;
		updateTableSize(DEFAULT_TABLE_SIZE - 1);
		System.out.println("\tАнализ быстрого численного метода решения уравнения Лапласа");
		requestNewTableSize();
	}

	public static void requestNewApproximations(){
		System.out.print("\nВведите новое значение аппрокисмации: ");
		int tmp = scanner.nextInt();
		if (tmp < 1)
			System.out.println("Аппроксимация не может быть меньше 1");
		else
			approximations = tmp;
	}

	public static void requestNewTableSize(){
		System.out.print("\nВведите размерность исследуемой области: ");
		int tmp = scanner.nextInt() - 1;
		if (tmp < DEFAULT_TABLE_SIZE - 1)
			System.out.printf("Размерность не может быть меньше %d%n", DEFAULT_TABLE_SIZE);
		else
			updateTableSize(tmp);
	}

	public static void updateTableSize(int newTableSize){
		tableSize = newTableSize;
		updateTables();
	}

	public static void updateTables(){
		initArrays();
		computeValueStep();
		computeBoardValues();
		F = 0;
		setupTableOfBoundaryConditions();
		laplasCorrectSolution();
	}

	public static void initArrays(){
		initValuesTablesArrays();
		initBoundArrays();
	}
	public static void initValuesTablesArrays(){
		U = new double[tableSize + 1][tableSize + 1];
		UT = new double[tableSize + 1][tableSize + 1];
		UF1 = new double[tableSize + 1][tableSize + 1];
	}
	public static void initBoundArrays(){
		Yup = new double[tableSize + 1];
		Ydown = new double[tableSize + 1];
		Xleft = new double[tableSize + 1];
		Xright = new double[tableSize + 1];
	}

	public static void computeValueStep(){
		ΔhX = (double) (B - A) / tableSize;
		ΔhY = (double) (D - C) / tableSize;
	}
	public static void computeBoardValues(){
		for (int i = 0; i <= tableSize; i++) {
			// Лаплас 1
			Yup[i] = Math.cos(C + ΔhY * i);
			Ydown[i] = Math.exp(1) * Math.cos(C + ΔhY * i);
			Xleft[i] = Math.exp(A + ΔhX * i);
			Xright[i] = Math.exp(A + ΔhX * i) * Math.cos(1);
			// Лаплас 2
			// Yup[i] = Math.sin(Math.PI * Y[i] / D);
			// Ydown[i] = 0;
			// Xleft[i] = Math.sin(Math.PI * X[i] / B);
			// Xright[i] = 0;
		}
	}

	public static void setupTableOfBoundaryConditions(){
		U = MatrixOfBoundaryConditions(U);
		UT = MatrixOfBoundaryConditions(UT);
		UF1 = MatrixOfBoundaryConditions(UF1);
	}
	public static void laplasCorrectSolution(){
		// Лаплас 1
		for (int i = 1; i < tableSize; i++)
			for (int j = 1; j < tableSize; j++)
				UT[i][j] = Math.cos(A + ΔhX * j) * Math.exp(C + ΔhY * i);
		// Лаплас 2
//		for (int i = 1; i < tableSize; i++)
//			for (int j = 1; j < tableSize; j++)
//				UT[i][j] = Math.sinh(Math.PI * (B - A + ΔhX * j) / D) / Math.sinh(Math.PI * B / D)
//					* Math.sin(Math.PI * (C + ΔhY * i) / D)
//					+ Math.sinh(Math.PI * (D - (C + ΔhY * i)) / B)
//					/ Math.sinh(Math.PI * D / B) * Math.sin(Math.PI * (A + ΔhX * j) / B);
	}

	public static void main(String[] args) {
		/**
		 * tableSize - количество шагов, n1 - количество шагов для вычисления интеграла Фурье,
		 * N, T - количество уравнений в СЛАУ для уравнения теплопроводности,
		 * a - коэффициент пропорциональности, l - длина, h - высота,
		 * A, B - границы исследуемого отрезка, t - исследуемое время,
		 * Δh, Δt - шаги по пространственной и временной/температурной координатам,
		 * λ - универсальная переменная начальных условий,
		 * calculationError - средняя квадратическая ошибка,
		 * X - массив значений пространственных координат,
		 * ΔU - вспомогательный массив значений искомой функции для записи при заданном
		 * времени, а также для вычисления средней квадратической ошибки,
		 * UT - массив значений искомой функции, вычисленных аналитическим методом,
		 * U - массив значений искомой функции, вычисленных численным методом,
		 * M - массив для составления трёхдиагональной матрицы,
		 */
		Menu menu = createMenu();
		menu.runInfinityLoop();
		System.out.println("Программа завершена.");
		scanner.close();
	}
	public static Menu createMenu(){
		return new Menu().withExitOnZero()
				.withMenuItem(1, "Решение уравнения Лапласа быстрым численным методом.", FNSoDE1::fastMethod)
				.withMenuItem(2, "Решение уравнения Лапласа численным методом конечных разностей.", FNSoDE1::finiteDifferenceMethod)
				.withMenuItem(3, "Включить / выключить вывод таблиц.", FNSoDE1::changePrintTableMod)
				.withMenuItem(4, "Включить / выключить режим многопоточной обработки.", FNSoDE1::changeMultiThreadMod)
				.withMenuItem(5, "Изменить размерность.", FNSoDE1::requestNewTableSize)
				.withMenuItem(6, "Изменить аппроксимацию.", FNSoDE1::requestNewApproximations);
	}

	public static void fastMethod(){
		long millis = System.currentTimeMillis();
		calculationError = 0;
		U = MatrixOfBoundaryConditions(new double[tableSize + 1][tableSize + 1]);
		for (int i = 1; i <= tableSize; i++) {
			for (int j = 1; j < tableSize; j++) {
				if (i != tableSize)
					UF1[i][j] = (UF1[i - 1][j] + UF1[i - 1][j - 1] + UF1[i][j - 1]) / 3;
				if (i != 1)
					U[i - 1][j] += (UF1[i - 2][j] + UF1[i][j] + UF1[i - 1][j - 1] + UF1[i - 1][j + 1]) / 16;
			}
		}
		for (int i = 1; i <= tableSize; i++) {
			for (int j = tableSize - 1; j >= 1; j--) {
				if (i != tableSize)
					UF1[i][j] = (UF1[i - 1][j] + UF1[i - 1][j + 1] + UF1[i][j + 1]) / 3;
				if (i != 1)
					U[i - 1][j] += (UF1[i - 2][j] + UF1[i][j] + UF1[i - 1][j - 1] + UF1[i - 1][j + 1]) / 16;
			}
		}
		for (int i = tableSize - 1; i >= 0; i--) {
			for (int j = 1; j < tableSize; j++) {
				if (i != 0)
					UF1[i][j] = (UF1[i + 1][j] + UF1[i + 1][j - 1] + UF1[i][j - 1]) / 3;
				if (i != tableSize - 1)
					U[i + 1][j] += (UF1[i][j] + UF1[i + 2][j] + UF1[i + 1][j - 1] + UF1[i + 1][j + 1]) / 16;
			}
		}
		for (int i = tableSize - 1; i >= 0; i--) {
			for (int j = tableSize - 1; j >= 1; j--) {
				if (i != 0)
					UF1[i][j] = (UF1[i + 1][j] + UF1[i + 1][j + 1] + UF1[i][j + 1]) / 3;
				if (i != tableSize - 1)
					U[i + 1][j] += (UF1[i][j] + UF1[i + 2][j] + UF1[i + 1][j - 1] + UF1[i + 1][j + 1]) / 16;
			}
		}
		System.out.printf("4 начальных приближения выполены: %d%n", System.currentTimeMillis() - millis);
		for (int k = 1, nCells = 50; k <= approximations; k++) {
			progressbarShow(k, approximations, nCells, "Аппроксимаций выполнено ");
			for (int i = 1; i < tableSize; i++) {
				for (int j = 1; j < tableSize; j++) {
					if (k % 2 != 0) {
						UF1[i][j] = (U[i - 1][j] + U[i + 1][j] + U[i][j - 1] + U[i][j + 1]) / 4;
					} else {
						U[i][j] = (UF1[i - 1][j] + UF1[i + 1][j] + UF1[i][j - 1] + UF1[i][j + 1]) / 4;
					}
				}
			}
		}
		if (approximations % 2 != 0) {
			U = UF1;
		}
		System.out.println("Решение уравнения Лапласа быстрым численным методом:");
		PrintTable(U);
		System.out.println("Точное решение уравнения Лапласа:");
		PrintTable(UT);
		for (int i = 1; i < tableSize; i++) {
			for (int j = 1; j < tableSize; j++) {
				calculationError += Math.abs((UT[i][j] - U[i][j]) / UT[i][j] * 100);
			}
		}
		calculationError = calculationError / Math.pow((tableSize - 1), 2);
		System.out.printf("Погрешность вычислений: %.6f%%%n", calculationError);
		System.out
				.println("Время проведённых вычислений: " + ((System.currentTimeMillis() - millis)) + " мс.");
	}

	public static void finiteDifferenceMethod(){
		long millis = System.currentTimeMillis();
		calculationError = 0;
		FiniteDifferenceMethodLP();
		System.out.println("Решение уравнения Лапласа численным методом конечных разностей:");
		PrintTable(U);
		System.out.println("Точное решение уравнения Лапласа:");
		PrintTable(UT);
		for (int i = 1; i < tableSize; i++) {
			for (int j = 1; j < tableSize; j++) {
				calculationError += Math.abs((UT[i][j] - U[i][j]) / UT[i][j] * 100);
			}
		}
		calculationError = calculationError / Math.pow((tableSize - 1), 2);
		System.out.printf("Погрешность вычислений в процентах: %.6f%n", calculationError);
		System.out.printf("Время проведённых вычислений: %d мс.%n", ((System.currentTimeMillis() - millis)));
	}

	public static void FiniteDifferenceMethodLP() { // метод конечных разностей
		int p = (tableSize - 1) * (tableSize - 1);
		double[][] M = new double[p][p + 1];
		double[] ΔU = new double[p];
		for (int i = 0; i < tableSize - 1; i++) {
			for (int j = i * (tableSize - 1); j < (i + 1) * (tableSize - 1); j++) {
				M[j][j] = 4;
				if (j % (tableSize - 1) != tableSize - 2) {
					M[j + 1][j] = -1;
					M[j][j + 1] = -1;
				}
				if (i != (tableSize - 2)) {
					M[j][j + tableSize - 1] = -1;
				}
				if (i != 0) {
					M[j][j - tableSize + 1] = -1;
				}
			}
		}
		M[0][p] = Xleft[1] + Yup[1];
		M[tableSize - 2][p] = Xright[1] + Yup[tableSize - 1];
		M[p - tableSize + 1][p] = Xleft[tableSize - 1] + Ydown[1];
		M[p - 1][p] = Xright[tableSize - 1] + Ydown[tableSize - 1];
		for (int i = 1; i < tableSize - 2; i++) {
			M[i][p] = Yup[i + 1];
			M[i * (tableSize - 1)][p] = Xleft[i + 1];
			M[i * (tableSize - 1) + (tableSize - 2)][p] = Xright[i + 1];
			M[i + p - tableSize + 1][p] = Ydown[i + 1];
		}
		if (F != 0) {
			for (int i = 0; i < p; i++) {
				M[i][p] -= ΔhX * ΔhY * F;
			}
		}
		printTableG(M,p);
		methodGauss(M, p, ΔU);
		for (int i = 1; i < tableSize; i++) {
			System.arraycopy(ΔU, (i - 1) * (tableSize - 1), U[i], 1, tableSize - 1);
		}
	}

	public static double[][] MatrixOfBoundaryConditions(double[][] matrix) { // присвоение таблицам граничных значений
		for (int i = 0; i <= tableSize; i++) {
			matrix[C][i] = Yup[i];
			matrix[D * tableSize][i] = Ydown[i];
			matrix[i][A] = Xleft[i];
			matrix[i][B * tableSize] = Xright[i];
		}
		return matrix;
	}

	public static void changePrintTableMod(){
		printingTable = !printingTable;
		System.out.println("Текущий режим вывода таблиц: " + printingTable);
	}
	public static void PrintTable(double[][] matrix) { // функция для вывода таблицы значений
		if (!printingTable)
			return;
		StringBuilder tableView = new StringBuilder();
		tableView.append("╔═══════╦");
		for (int i = 0; i < tableSize; i++) {
			tableView.append("══════════╦");
		}
		tableView.append("══════════╗\n║ U(x,y)");
		for (int i = 0; i < tableSize + 1; i++) { // цикл, заполняющий шапку таблицы
			tableView.append(String.format("║   %.2f   ", A + ΔhX * i));
		}
		tableView.append("║\n╠═══════╬");
		for (int i = 0; i < tableSize; i++) {
			tableView.append("══════════╬");
		}
		tableView.append("══════════╣\n");
		for (int i = 0; i < tableSize + 1; i++) { // цикл, проходящий строки таблицы
			tableView.append(String.format("║ %.2f\t", C + ΔhY * i));
			for (int j = 0; j < tableSize + 1; j++) { // цикл, заполняющий ячейки текущей строки
				tableView.append(String.format("║ %.6f ", matrix[i][j]));
			}
			tableView.append("║\n");
		}
		tableView.append("╚═══════╩");
		for (int i = 0; i < tableSize; i++) {
			tableView.append("══════════╩");
		}
		tableView.append("══════════╝\n");
		System.out.println(tableView);
	}

	public static void methodGauss(double[][] A, int tableSize, double[] X) { // функция, реализующая метод Гаусса
		double[] B = new double[tableSize];
		for (int i = 0; i < tableSize; i++) {
			B[i] = A[i][tableSize];
		}
		// прямой ход (или приведение матрицы к треугольному виду)

		for (int p = 0, nCells = 50; p < tableSize; p++) { // с помощью цикла перебираем все уравнения
			progressbarShow(p, tableSize, nCells, "Обработано строк ");
			int max = p; // задаём максимальное значение первого ненулевого коэффициента в уравнении
			for (int i = p + 1; i < tableSize; i++) { // пробегаем по столбцам и ищем максимальное значение
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
			for (int i = p + 1; i < tableSize; i++) { // цикл, преобразующий строку так, что её первый элемент, считая над
															// первым
															// ненулевым элементом, становится равен нулю
				double alpha = A[i][p] / A[p][p]; // получение коэффициента, на который будет умножаться строка,
																// отнимаемая
																// от всех последующих
				B[i] -= alpha * B[p]; // отдельно производим вычитание для свободных членов
				for (int j = p; j < tableSize; j++) {
					A[i][j] -= alpha * A[p][j]; // производим вычитание коэффициента, умноженного на строку так, что её
															// первый элемент, считая над первым ненулевым элементом, становится
															// равен
															// нулю
				}
			}
		}
		// обратный ход
		for (int i = tableSize - 1; i >= 0; i--) { // обратный ход осуществим через два цикла; пройдём все строки СЛАУ в обратном
														// порядке
			double sum = 0.0; // вводим счётчик
			for (int j = i + 1; j < tableSize; j++) { // цикл заключается в получении в счётчик значения, при котором при
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
	public static void printTableG(double A[][], int n) {
		StringBuilder tableView = new StringBuilder();
		for (int i = 0; i < n; i++) { // цикл, проходящий строки матрицы
			for (int j = 0; j < n; j++) { // цикл, заполняющий ячейки текущей строки
				tableView.append(String.format("║  %.3f\t", A[i][j]));
			}
			tableView.append(String.format("║  %.3f\t", A[i][n]));
		}
		System.out.println(tableView);
	}

	public static void progressbarShow(int nTaskCompleted, int nTaskCount, int nCells, String comment){
		StringBuilder stringBuilder = new StringBuilder();
		stringBuilder.append(comment)
				.append(" |");
		int progress = nCells * nTaskCompleted / nTaskCount ;
		for (int i = 0; i < progress; i++)
			stringBuilder.append("=");
		for (int i = progress; i < nCells; i++)
			stringBuilder.append("Х");
		stringBuilder.append("| ")
				.append(nTaskCompleted)
				.append(" / ")
				.append(nTaskCount)
				.append("\r");
		System.out.print(stringBuilder);
	}
}