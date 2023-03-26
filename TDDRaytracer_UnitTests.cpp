#include "pch.h"
#include "CppUnitTest.h"
#include "C:\Users\strai\source\TDD_raytracer\TDD_Raytracer\ArithmeticStructures.h"
#include "C:\Users\strai\source\TDD_raytracer\TDD_Raytracer\GeometricStructures.h"
#include "C:\Users\strai\source\TDD_raytracer\TDD_Raytracer\SceneObject.h"
#include "C:\Users\strai\source\TDD_raytracer\TDD_Raytracer\Ray.h"
#include "C:\Users\strai\source\TDD_raytracer\TDD_Raytracer\Canvas.h"
#include "C:\Users\strai\source\TDD_raytracer\TDD_Raytracer\PPMWriter.h"
#include <fstream>
#include <filesystem>
#include <array>
#include <cmath>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace TDDRaytracerUnitTests
{
	TEST_CLASS(TDDRaytracerUnitTests)
	{
	public:
		
		TEST_METHOD(ArithmeticStructures_IsVectorTest)
		{
			ArithmeticStructures aS;
			aS.setVector(1.3, 2.3, -4.0);
			
			constexpr float expectedWPoint{ 0.0 };
			const auto &[x,y,z,w] = aS.getVector();

			Assert::AreEqual(expectedWPoint, w, 0.0001f, L"Vectors homogeneous coordinate must be 0.0");
		}

		TEST_METHOD(ArithmeticStructures_IsPointTest)
		{
			ArithmeticStructures aS;
			aS.setPoint(1.3, 2.3, -4.0);

			constexpr float expectedWVector{ 1.0 };
			const auto& [x, y, z, w] = aS.getPoint();

			Assert::AreEqual(expectedWVector, w, 0.0001f, L"Points homogeneous coordinate must be 1.0");
		}

		TEST_METHOD(ArithmeticStructures_EqualityTest)
		{
			ArithmeticStructures aS1;
			aS1.setPoint(1.3, 2.3, -4.0);
			aS1.setVector(1.3, 2.3, -4.0);
			ArithmeticStructures aS2;
			aS2.setPoint(1.3, 2.3, -4.0);
			aS2.setVector(1.3, 2.3, -4.0);

			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(aS1.getPoint(), aS2.getPoint()));
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(aS1.getVector(), aS2.getVector()));

			Assert::IsFalse(ArithmeticStructures::coordinatesAreEqual(aS1.getVector(), aS2.getPoint()));
			
		}

		TEST_METHOD(ArithmeticStructures_AddCoordinatesTest)
		{
			ArithmeticStructures aS1;
			aS1.setPoint(1.0, 1.0, 1.0);
			aS1.setVector(1.0, 1.0, 1.0);
			ArithmeticStructures aS2;
			aS2.setPoint(1.0, 1.0, 1.0);
			aS2.setVector(1.0, 1.0, 1.0);
			ArithmeticStructures expectedAS;
			expectedAS.setPoint(2.0, 2.0, 2.0);
			expectedAS.setVector(2.0, 2.0, 2.0);

			// point + vector results in point (w == 1)
			ArithmeticStructures::HomogenousCoordinates sumOfCoordinatesIsPoint{ ArithmeticStructures::addCoordinates(aS1.getPoint(), aS2.getVector()) };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(sumOfCoordinatesIsPoint, expectedAS.getPoint()));

			// vector + vector results in vector (w == 0)
			ArithmeticStructures::HomogenousCoordinates sumOfCoordinatesIsVector{ ArithmeticStructures::addCoordinates(aS1.getVector(), aS2.getVector()) };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(sumOfCoordinatesIsVector, expectedAS.getVector()));

		}

		TEST_METHOD(ArithmeticStructures_SubtractCoordinatesTest)
		{
			ArithmeticStructures::HomogenousCoordinates point_1{ 3.0,2.0,1.0,1.0 };
			ArithmeticStructures::HomogenousCoordinates point_2{ 5.0,6.0,7.0,1.0 };
			const ArithmeticStructures::HomogenousCoordinates expectedVector_afterSubtraction{ -2.0,-4.0,-6.0,0.0 };
			// subtracting point from point results in expected vector
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedVector_afterSubtraction, ArithmeticStructures::subtractCoordinates(point_1, point_2)));

			// "transform" the point into a vector
			auto& [x_point2, y_point2, z_point2, w_point2] = point_2;
			w_point2 = 0.0;
			// subtracting point from vector results in expected vector
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedVector_afterSubtraction, ArithmeticStructures::subtractCoordinates(point_1, point_2)));

			// "transform" the point into a vector
			auto& [x_point1, y_point1, z_point1, w_point1] = point_1;
			w_point1 = 0.0;
			// subtracting vector from vector results in expected vector
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedVector_afterSubtraction, ArithmeticStructures::subtractCoordinates(point_1, point_2)));
		}

		TEST_METHOD(ArithmeticStructures_MultiplyScalarTest)
		{
			constexpr float s{ 2.0 };
			ArithmeticStructures aS;
			aS.setVector( 1.0,0.0,0.0 );
			aS.setPoint(0.2, 0.3, 0.4);
			ArithmeticStructures expectedValue;
			expectedValue.setVector(2.0, 0.0, 0.0);
			expectedValue.setPoint(0.4, 0.6, 0.8);
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(ArithmeticStructures::multiplyWithScalar(aS.getVector(), s), expectedValue.getVector()));
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(ArithmeticStructures::multiplyWithScalar(aS.getPoint(), s), expectedValue.getPoint()));
		}

		TEST_METHOD(ArithmeticStructures_MultiplyVectorTest)
		{
			ArithmeticStructures aS1;
			aS1.setVector(1.0, 0.2, 0.4);
			ArithmeticStructures aS2;
			aS2.setVector(0.9, 1.0, 0.1);
			ArithmeticStructures expectedValue;
			expectedValue.setVector(0.9, 0.2, 0.04);
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(ArithmeticStructures::multiplyWithVector(aS1.getVector(), aS2.getVector()), expectedValue.getVector()));
		}

		TEST_METHOD(ArithmeticStructure_MagnitudeTest)
		{
			ArithmeticStructures aS;
			aS.setPoint(1.3, 2.3, -4.0);
			aS.setVector(1.0, 0.0, 0.0);

			constexpr float epsilon{ 0.0001 };
			float expectedMagnitude{ 1.0 };

			Assert::AreEqual(expectedMagnitude, aS.magnitude(), epsilon);
			aS.setVector(0.0, 1.0, 0.0);
			Assert::AreEqual(expectedMagnitude, aS.magnitude(), epsilon);

			expectedMagnitude = sqrtf(14);
			aS.setVector(1.0, 2.0, 3.0);
			Assert::AreEqual(expectedMagnitude, aS.magnitude(), epsilon);

			aS.setVector(-1.0, 2.0,-3.0);
			Assert::AreEqual(expectedMagnitude, aS.magnitude(), epsilon);

		}

		TEST_METHOD(ArithmeticStructure_NormalizationTest)
		{
			ArithmeticStructures aS;
			aS.setVector(1.3, 2.3, -4.0);

			constexpr float epsilon{ 0.0001 };
			constexpr float expectedMagnitude{ 1.0 };

			const auto [x, y, z, w] = aS.getNormalizedVector();

			 ArithmeticStructures n;
			 n.setVector(x, y, z);

			Assert::AreEqual(expectedMagnitude, n.magnitude(), epsilon);

		}

		TEST_METHOD(ArithmeticStructure_DotProductTest)
		{
			ArithmeticStructures aS1;
			aS1.setVector(1.0, 2.0, 3.0);

			ArithmeticStructures aS2;
			aS2.setVector(2.0, 3.0, 4.0);

			constexpr float epsilon{ 0.0001 };
			constexpr float expectedValue{ 20.0 };


			Assert::AreEqual(expectedValue, ArithmeticStructures::dotProduct(aS1.getVector(), aS2.getVector()), epsilon);

		}

		TEST_METHOD(ArithmeticStructure_CrossProductTest)
		{
			ArithmeticStructures aS1;
			aS1.setVector(1.0, 2.0, 3.0);

			ArithmeticStructures aS2;
			aS2.setVector(2.0, 3.0, 4.0);

			ArithmeticStructures expectedValue;
			expectedValue.setVector(-1.0, 2.0, -1.0);

			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedValue.getVector(), ArithmeticStructures::crossProduct(aS1.getVector(), aS2.getVector())));

			expectedValue.setVector(1.0, -2.0, 1.0);
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedValue.getVector(), ArithmeticStructures::crossProduct(aS2.getVector(), aS1.getVector())));
		}

		TEST_METHOD(ArithmeticStructure_MatrixTypeTest)
		{
			ArithmeticStructures::row2x2 m0_2x2{ {0.0, 0.0} }, m1_2x2{ {0.0, 0.0} };
			ArithmeticStructures::Matrix2x2 m22{m0_2x2, m1_2x2};
			Assert::IsTrue(ArithmeticStructures::MatrixType::Matrix2x2 == m22.getType());

			ArithmeticStructures::row3x3 m0_3x3{ {0.0, 0.0, 0.0} }, m1_3x3{ {0.0, 0.0, 0.0} }, m2_3x3{ {0.0,0.0,0.0} };
			ArithmeticStructures::Matrix3x3 m33{ m0_3x3, m1_3x3, m2_3x3 };
			Assert::IsTrue(ArithmeticStructures::MatrixType::Matrix3x3 == m33.getType());

			ArithmeticStructures::row4x4 m0_4x4{ {0.0,0.0, 0.0, 0.0} }, m1_4x4{ {0.0, 0.0, 0.0, 0.0} }, m2_4x4{ {0.0, 0.0,0.0,0.0} }, m3_4x4{ {0.0,0.0,0.0,0.0} };
			ArithmeticStructures::Matrix4x4 m44{m0_4x4, m1_4x4, m2_4x4, m3_4x4};
			Assert::IsTrue(ArithmeticStructures::MatrixType::Matrix4x4 == m44.getType());
		}

		TEST_METHOD(ArithmeticStructure_MatrixSizeTest)
		{
			ArithmeticStructures::row2x2 m0_2x2{ {0.0, 0.0} }, m1_2x2{ {0.0, 0.0} };
			ArithmeticStructures::Matrix2x2 m22{ m0_2x2, m1_2x2 };
			Assert::AreEqual(4, m22.getMatrixSize());

			ArithmeticStructures::row3x3 m0_3x3{ {0.0, 0.0, 0.0} }, m1_3x3{ {0.0, 0.0, 0.0} }, m2_3x3{ {0.0,0.0,0.0} };
			ArithmeticStructures::Matrix3x3 m33{ m0_3x3, m1_3x3, m2_3x3 };
			Assert::AreEqual(9, m33.getMatrixSize());

			ArithmeticStructures::row4x4 m0_4x4{ {0.0, 0.0, 0.0, 0.0} }, m1_4x4{ {0.0, 0.0, 0.0, 0.0} }, m2_4x4{ {0.0,0.0,0.0, 0.0} }, m3_4x4{ {0.0,0.0,0.0,0.0} };
			ArithmeticStructures::Matrix4x4 m44{ m0_4x4, m1_4x4, m2_4x4, m3_4x4 };
			Assert::AreEqual(16, m44.getMatrixSize());
		}

		TEST_METHOD(ArithmeticStructure_MatrixComparisonTest)
		{
			ArithmeticStructures::row2x2 m0_2x2{ {0.0, 1.0} }, m1_2x2{ {0.0, 0.0} };
			ArithmeticStructures::Matrix2x2 matrix1_2x2{ m0_2x2, m1_2x2 }, matrix2_2x2{ m0_2x2, m1_2x2 };
			// detect when matrices are equal
			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_2x2(matrix1_2x2, matrix2_2x2));
			// detect when matrices are not equal
			matrix2_2x2.setMatrixData(1, ArithmeticStructures::row2x2{ 3.0,0.0 });
			Assert::IsFalse(ArithmeticStructures::matricesAreEqual_2x2(matrix1_2x2, matrix2_2x2));

			ArithmeticStructures::row3x3 m0_3x3{ { 0.0, 3.0, 0.0} }, m1_3x3{ { 0.0, 0.0, 0.0} }, m2_3x3{ {0.0,0.0, 0.0}  };
			ArithmeticStructures::Matrix3x3 matrix1_3x3{ m0_3x3, m1_3x3, m2_3x3}, matrix2_3x3{ m0_3x3, m1_3x3, m2_3x3 };
			// detect when matrices are equal
			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_3x3(matrix1_3x3, matrix2_3x3));
			// detect when matrices are not equal
			matrix2_3x3.setMatrixData(1, ArithmeticStructures::row3x3{ 3.0,0.0,0.0 });
			Assert::IsFalse(ArithmeticStructures::matricesAreEqual_3x3(matrix1_3x3, matrix2_3x3));

			ArithmeticStructures::row4x4 m0_4x4{ {0.0, 0.0, 0.0, 0.0} }, m1_4x4{ {0.0, 0.0, 0.0, 0.0} }, m2_4x4{ {0.0,0.0,0.0, 0.0} }, m3_4x4{ {0.0,0.0,0.0,0.0} };
			ArithmeticStructures::Matrix4x4 matrix1_4x4{ m0_4x4, m1_4x4, m2_4x4, m3_4x4 }, matrix2_4x4{ m0_4x4, m1_4x4, m2_4x4, m3_4x4 };
			// detect when matrices are equal
			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(matrix1_4x4, matrix2_4x4));
			// detect when matrices are not equal
			matrix2_4x4.setMatrixData(1, ArithmeticStructures::row4x4{ 3.0,0.0,0.0,0.0 });
			Assert::IsFalse(ArithmeticStructures::matricesAreEqual_4x4(matrix1_4x4, matrix2_4x4));

		}

		TEST_METHOD(ArithmeticStructure_MatrixValuesTest)
		{
			constexpr float expectedValue{ 1.0 };
			ArithmeticStructures::row2x2 m0_2x2{ {0.0, expectedValue} }, m1_2x2{ {0.0, 0.0} };
			ArithmeticStructures::Matrix2x2 m2x2{ m0_2x2, m1_2x2 };
			Assert::AreEqual(expectedValue, m2x2.getElement(0, 1));


			ArithmeticStructures::row3x3 m0_3x3{ { 0.0, 3.0, 0.0} }, m1_3x3{ { 0.0, expectedValue, 0.0} }, m2_3x3{ {0.0,0.0, 0.0} };
			ArithmeticStructures::Matrix3x3 m3x3{ m0_3x3, m1_3x3, m2_3x3 };
			Assert::AreEqual(expectedValue, m3x3.getElement(1, 1));

			ArithmeticStructures::row4x4 m0_4x4{ {0.0, 0.0, 0.0, 0.0} }, m1_4x4{ {0.0, 0.0, 0.0, 0.0} }, m2_4x4{ {0.0,0.0,0.0, 0.0} }, m3_4x4{ {0.0,0.0,expectedValue,0.0} };
			ArithmeticStructures::Matrix4x4 m4x4{ m0_4x4, m1_4x4, m2_4x4, m3_4x4 };
			Assert::AreEqual(expectedValue, m4x4.getElement(3, 2));
		}

		TEST_METHOD(ArithmeticStructure_MatrixMatrixMultiplicationTest)
		{
			// test matrix X matrix multiplication
			// values for matrix1
			ArithmeticStructures::row4x4 m0_4x4_1{ {1.0, 2.0, 3.0, 4.0} }, m1_4x4_1{ {5.0, 6.0, 7.0, 8.0} }, m2_4x4_1{ {9.0,8.0,7.0, 6.0} }, m3_4x4_1{ {5.0,4.0,3.0,2.0} };
			// values for matrix2
			ArithmeticStructures::row4x4 m0_4x4_2{ {-2.0, 1.0, 2.0, 3.0} }, m1_4x4_2{ {3.0, 2.0, 1.0, -1.0} }, m2_4x4_2{ {4.0,3.0,6.0, 5.0} }, m3_4x4_2{ {1.0,2.0,7.0,8.0} };
			// values for expected result
			ArithmeticStructures::row4x4 m0_expected{ {20.0, 22.0, 50.0, 48.0} }, m1_expected{ {44.0, 54.0, 114.0, 108.0} }, m2_expected{ {40.0,58.0,110.0, 102.0} }, m3_expected{ {16.0,26.0,46.0,42.0} };
			
			ArithmeticStructures::Matrix4x4 
				m1{ m0_4x4_1, m1_4x4_1, m2_4x4_1, m3_4x4_1 }, 
				m2{ m0_4x4_2, m1_4x4_2, m2_4x4_2, m3_4x4_2 }, 
				expectedResult{ m0_expected, m1_expected, m2_expected, m3_expected };
			
			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(expectedResult, ArithmeticStructures::multiplyMatrices(m1,m2)));

			
		}

		TEST_METHOD(ArithmeticStructure_MatrixIdentityTest)
		{
			ArithmeticStructures::row4x4 m0{ {1.0, 2.0, 3.0, 4.0} }, m1{ {2.0, 4.0, 4.0, 2.0} }, m2{ {8.0,6.0,4.0, 1.0} }, m3{ {0.0,0.0,0.0,1.0} };
			ArithmeticStructures::Matrix4x4 m{ m0, m1, m2, m3 };
			ArithmeticStructures::Matrix4x4 expectedResult{ m.getRowM(0),m.getRowM(1),m.getRowM(2),m.getRowM(3) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(expectedResult, ArithmeticStructures::multiplyMatrices(m, ArithmeticStructures::getIdentityMatrix())));


		}

		TEST_METHOD(ArithmeticStructure_MatrixTranspositionTest)
		{
			ArithmeticStructures::row4x4 m0{ {0.0, 9.0, 3.0, 0.0} }, m1{ {9.0, 8.0, 0.0, 8.0} }, m2{ {1.0,8.0,5.0, 3.0} }, m3{ {0.0,0.0,5.0,8.0} };
			ArithmeticStructures::Matrix4x4 m{ m0, m1, m2, m3 };
			ArithmeticStructures::row4x4 expectedResult_m0{ {0.0, 9.0, 1.0, 0.0} }, expectedResult_m1{ {9.0, 8.0, 8.0, 0.0} }, expectedResult_m2{ {3.0,0.0,5.0, 5.0} }, expectedResult_m3{ {0.0,8.0,3.0,8.0} };
			ArithmeticStructures::Matrix4x4 expectedResult{ expectedResult_m0 , expectedResult_m1,expectedResult_m2,expectedResult_m3 };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(expectedResult, ArithmeticStructures::transposeMatrix(m)));

			ArithmeticStructures::Matrix4x4 identiyt{ 
				ArithmeticStructures::getIdentityMatrix().getRowM(0) , 
				ArithmeticStructures::getIdentityMatrix().getRowM(1),
				ArithmeticStructures::getIdentityMatrix().getRowM(2),
				ArithmeticStructures::getIdentityMatrix().getRowM(3) };
			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(ArithmeticStructures::getIdentityMatrix(), ArithmeticStructures::transposeMatrix(identiyt)));


		}

		TEST_METHOD(ArithmeticStructure_MatrixCalculateDeterminantTest)
		{
			ArithmeticStructures::row2x2 m0{ {1.0, 5.0} }, m1{ {-3.0, 2.0} };
			ArithmeticStructures::Matrix2x2 m{ m0, m1};

			float expectedValue{ 17 };

			Assert::AreEqual(expectedValue, m.getDeterminant());

			ArithmeticStructures::row3x3 m0_3x3{ {1.0, 2.0, 6.0} }, m1_3x3{ {-5.0, 8.0, -4.0} }, m2_3x3{ {2.0,6.0,4.0} };
			ArithmeticStructures::Matrix3x3 m_3x3{ m0_3x3 , m1_3x3,m2_3x3 };
			expectedValue = -196;

			Assert::AreEqual(expectedValue, m_3x3.getDeterminant());

			ArithmeticStructures::row4x4 m0_4x4{ {-2.0, -8.0, 3.0, 5.0} }, m1_4x4{ {-3.0, 1.0, 7.0, 3.0} }, m2_4x4{ {1.0,2.0,-9.0, 6.0} }, m3_4x4{ {-6.0,7.0,7.0,-9.0} };
			ArithmeticStructures::Matrix4x4 m_4x4{ m0_4x4, m1_4x4, m2_4x4, m3_4x4 };
			expectedValue = -4071;

			Assert::AreEqual(expectedValue, m_4x4.getDeterminant());

			
		}

		TEST_METHOD(ArithmeticStructure_MatrixIsInvertibleTest)
		{
			// this matrix is expected to be invertible
			ArithmeticStructures::row4x4 m0_1{ {6.0, 4.0, 4.0, 4.0} }, m1_1{ {5.0, 5.0, 7.0, 6.0} }, m2_1{ {4.0,-9.0,3.0, -7.0} }, m3_1{ {9.0,1.0,7.0,-6.0} };
			ArithmeticStructures::Matrix4x4 m_1{ m0_1, m1_1, m2_1, m3_1 };
			Assert::IsTrue(m_1.isInvertible());

			// this matrix is not expected to be invertible
			ArithmeticStructures::row4x4 m0_2{ {-4.0, 2.0, -2.0, -3.0} }, m1_2{ {9.0, 6.0, 2.0, 6.0} }, m2_2{ {0.0,-5.0,1.0, -5.0} }, m3_2{ {0.0,0.0,0.0,0.0} };
			ArithmeticStructures::Matrix4x4 m_2{ m0_2, m1_2, m2_2, m3_2 };
			Assert::IsFalse(m_2.isInvertible());
		}

		TEST_METHOD(ArithmeticStructure_MatrixInverseTest)
		{
			// this matrix is expected to be invertible
			ArithmeticStructures::row4x4 m0_1{ {-5.0, 2.0, 6.0, -8.0} }, m1_1{ {1.0, -5.0, 1.0, 8.0} }, m2_1{ {7.0,7.0,-6.0, -7.0} }, m3_1{ {1.0,-3.0,7.0,4.0} };
			ArithmeticStructures::Matrix4x4 m{ m0_1, m1_1, m2_1, m3_1 };

			ArithmeticStructures::row4x4 
				m0_expected{ {0.21805, 0.45113, 0.24060, -0.04511} }, 
				m1_expected{ {-0.80827, -1.45677, -0.44361, 0.52068} },
				m2_expected{ {-0.07895,-0.22368,-0.05263, 0.19737} }, 
				m3_expected{ {-0.52256,-0.81397,-0.30075,0.30639} };
			ArithmeticStructures::Matrix4x4 expectedResult{ m0_expected, m1_expected, m2_expected, m3_expected };

			ArithmeticStructures::inverseMatrix(m);

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(expectedResult, m));


			// this matrix is not expected to be invertible, inversion doesnt change that matrix though
			ArithmeticStructures::row4x4 m0_2{ {-4.0, 2.0, -2.0, -3.0} }, m1_2{ {9.0, 6.0, 2.0, 6.0} }, m2_2{ {0.0,-5.0,1.0, -5.0} }, m3_2{ {0.0,0.0,0.0,0.0} };
			ArithmeticStructures::Matrix4x4 m_2{ m0_2, m1_2, m2_2, m3_2 };

			ArithmeticStructures::inverseMatrix(m_2);

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_2, m_2));

			// test whether multiplying with the inverse behaves as expected (i.e. A*B = C -> C* B^-1 = A
			ArithmeticStructures::row4x4 m0_A{ {3.0, -9.0, 7.0, 3.0} }, m1_A{ {3.0, -8.0, 2.0, -9.0} }, m2_A{ {-4.0,4.0,4.0, 1.0} }, m3_A{ {-6.0,5.0,-1.0,1.0} };
			ArithmeticStructures::Matrix4x4 m_A{ m0_A, m1_A, m2_A, m3_A };
			ArithmeticStructures::row4x4 m0_B{ {8.0, 2.0, 2.0, 2.0} }, m1_B{ {3.0, -1.0, 7.0, 0.0} }, m2_B{ {7.0,0.0,5.0, 4.0} }, m3_B{ {6.0,-2.0,0.0,5.0} };
			ArithmeticStructures::Matrix4x4 m_B{ m0_B, m1_B, m2_B, m3_B };

			auto m_C{ ArithmeticStructures::multiplyMatrices(m_A,m_B) };
			ArithmeticStructures::inverseMatrix(m_B);
			auto multiplicationWithInverseResult{ ArithmeticStructures::multiplyMatrices(m_C,m_B) };
			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_A, multiplicationWithInverseResult));
		}

		TEST_METHOD(ArithmeticStructure_MatrixCofactorMatrixTest)
		{
			// this matrix is expected to be invertible
			ArithmeticStructures::row4x4 m0_1{ {-5.0, 2.0, 6.0, -8.0} }, m1_1{ {1.0, -5.0, 1.0, 8.0} }, m2_1{ {7.0,7.0,-6.0, -7.0} }, m3_1{ {1.0,-3.0,7.0,4.0} };
			ArithmeticStructures::Matrix4x4 m{ m0_1, m1_1, m2_1, m3_1 };

			ArithmeticStructures::row4x4
				m0_2{ {116, -430, -42, -278} },
				m1_2{ {240, -775, -119, -433} },
				m2_2{ {128,-236,-28, -160} },
				m3_2{ {-24,277,105,163} };
			ArithmeticStructures::Matrix4x4 expectedResult{ m0_2, m1_2, m2_2, m3_2 };

			

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(expectedResult, ArithmeticStructures::getCofactorMatrix(m)));
		}


		TEST_METHOD(ArithmeticStructure_MatrixSubmatrixTest)
		{
			ArithmeticStructures::row4x4 m0{ {-6.0, 1.0, 1.0, 6.0} }, m1{ {-8.0, 5.0, 8.0, 6.0} }, m2{ {-1.0,0.0,8.0, 2.0} }, m3{ {-7.0,1.0,-1.0,1.0} };
			ArithmeticStructures::Matrix4x4 m{ m0, m1, m2, m3 };
			ArithmeticStructures::row3x3 expectedResult_m0{ {-6.0, 1.0, 6.0} }, expectedResult_m1{ {-8.0, 8.0, 6.0} }, expectedResult_m2{ {-7.0,-1.0,1.0} };
			ArithmeticStructures::Matrix3x3 expectedResult{ expectedResult_m0 , expectedResult_m1,expectedResult_m2};

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_3x3(expectedResult, ArithmeticStructures::getSubmatrixOf4x4Matrix(2,1,m)));



			ArithmeticStructures::row3x3 m0_3x3{ {1.0, 5.0, 0.0} }, m1_3x3{ {-3.0, 2.0, 7.0} }, m2_3x3{ {0.0,6.0,-3.0} };
			ArithmeticStructures::Matrix3x3 m_3x3{ m0_3x3, m1_3x3, m2_3x3};
			ArithmeticStructures::row2x2 expectedResult_m0_3x3{ {-3.0, 2.0} }, expectedResult_m1_3x3{ {0.0, 6.0} };
			ArithmeticStructures::Matrix2x2 expectedResult_3x3{ expectedResult_m0_3x3 , expectedResult_m1_3x3};

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_2x2(expectedResult_3x3, ArithmeticStructures::getSubmatrixOf3x3Matrix(0, 2, m_3x3)));
		}

		TEST_METHOD(ArithmeticStructure_MatrixMinorTest)
		{

			ArithmeticStructures::row3x3 m0_3x3{ {3.0, 5.0, 0.0} }, m1_3x3{ {2.0, -1.0, -7.0} }, m2_3x3{ {6.0,-1.0,5.0} };
			ArithmeticStructures::Matrix3x3 m_3x3{ m0_3x3, m1_3x3, m2_3x3 };

			constexpr float expectedResult{ 25.0 };
			constexpr int subMatrixRow{ 1 }, subMatrixColumn{ 0 };
			constexpr float epsilon{ 0.0001 };

			Assert::AreEqual(expectedResult, ArithmeticStructures::getSubmatrixOf3x3Matrix(subMatrixRow, subMatrixColumn, m_3x3).getDeterminant(), epsilon);

			Assert::AreEqual(  expectedResult, ArithmeticStructures::getMinor(subMatrixRow, subMatrixColumn, m_3x3), epsilon);
		}

		TEST_METHOD(ArithmeticStructure_MatrixCofactorTest)
		{

			ArithmeticStructures::row3x3 m0_3x3{ {3.0, 5.0, 0.0} }, m1_3x3{ {2.0, -1.0, -7.0} }, m2_3x3{ {6.0,-1.0,5.0} };
			ArithmeticStructures::Matrix3x3 m_3x3{ m0_3x3, m1_3x3, m2_3x3 };

			float expectedResult{ -25.0 };
			int subMatrixRow{ 1 }, subMatrixColumn{ 0 };
			constexpr float epsilon{ 0.0001 };

			Assert::AreEqual(expectedResult, ArithmeticStructures::getCofactor(subMatrixRow, subMatrixColumn, m_3x3), epsilon);
			subMatrixRow = 0; subMatrixColumn = 0;
			expectedResult = -12.0;
			Assert::AreEqual(expectedResult, ArithmeticStructures::getCofactor(subMatrixRow, subMatrixColumn, m_3x3), epsilon);
		}

		TEST_METHOD(ArithmeticStructure_MatrixTupleMultiplicationTest)
		{
			// test matrix X tuple multiplication
			ArithmeticStructures::row4x4 m0{ {1.0, 2.0, 3.0, 4.0} }, m1{ {2.0, 4.0, 4.0, 2.0} }, m2{ {8.0,6.0,4.0, 1.0} }, m3{ {0.0,0.0,0.0,1.0} };
			ArithmeticStructures::Matrix4x4 m{ m0, m1, m2, m3 };
			ArithmeticStructures::HomogenousCoordinates c{ 1.0,2.0,3.0,1.0 };
			ArithmeticStructures::HomogenousCoordinates expectedResult{18.0,24.0,33.0,1.0};
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult, ArithmeticStructures::multiplyMatrixWithTuple(m, c)));
		}

		TEST_METHOD(ArithmeticStructure_MatrixTranslationTest)
		{
			constexpr float shift_x{ 5.0 }, shift_y{ -3.0 }, shift_z{ 2 };

			// is the transformation matrix created as expected
			ArithmeticStructures::row4x4 m0_expected{ {1.0, 0.0, 0.0, shift_x} }, m1_expected{ {0.0, 1.0, 0.0, shift_y} }, m2_expected{ {0.0,0.0,1.0, shift_z} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };
			auto transformationMatrix{ ArithmeticStructures::getTranslationMatrix(shift_x,shift_y,shift_z) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));


			// does translation work as expected - M * p = p`
			ArithmeticStructures::HomogenousCoordinates originalPoint{ -3.0,4.0,5.0,1.0 };
			ArithmeticStructures::HomogenousCoordinates expectedResult_translatedPoint{ 2.0,1.0,7.0,1.0 };
			
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_translatedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));
			
			// and also test the inverse - M^-1 * p` = p
			ArithmeticStructures::inverseMatrix(transformationMatrix);
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(originalPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, expectedResult_translatedPoint)));
			
		}

		TEST_METHOD(ArithmeticStructure_MatrixScalingTest)
		{
			constexpr float scale_x{ 2.0 }, scale_y{ 3.0 }, scale_z{ 4.0 };

			// is the transformation matrix created as expected
			ArithmeticStructures::row4x4 m0_expected{ {scale_x, 0.0, 0.0,0.0} }, m1_expected{ {0.0, scale_y, 0.0, 0.0} }, m2_expected{ {0.0,0.0,scale_z, 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };
			auto transformationMatrix{ ArithmeticStructures::getScalingMatrix(scale_x,scale_y,scale_z) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does translation work as expected - M * p = p`
			ArithmeticStructures::HomogenousCoordinates originalPoint{ -4.0,6.0,8.0,1.0 };
			ArithmeticStructures::HomogenousCoordinates expectedResult_scaledPoint{ -8.0,18.0,32.0,1.0 };

			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_scaledPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));

			// also make sure that scaling of vectors works as expceted
			ArithmeticStructures::HomogenousCoordinates originalVector{ -4.0,6.0,8.0,0.0 };
			ArithmeticStructures::HomogenousCoordinates expectedResult_scaledVector{ -8.0,18.0,32.0,0.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_scaledVector, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalVector)));

			// and also test the inverse - M^-1 * p` = p
			ArithmeticStructures::inverseMatrix(transformationMatrix);
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(originalVector, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, expectedResult_scaledVector)));

			// also test reflection (aka "mirroring")
			constexpr float reflect_x{ -1.0 }, reflect_y{ 1.0 }, reflect_z{ 1.0 };
			ArithmeticStructures::HomogenousCoordinates originalPoint_2{ 2.0,3.0,4.0,1.0 };
			ArithmeticStructures::HomogenousCoordinates expectedResult_reflectedPoint{ -2.0,3.0,4.0,1.0 };
			auto transformationMatrix_reflection{ ArithmeticStructures::getScalingMatrix(reflect_x,reflect_y,reflect_z) };

			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_reflectedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix_reflection, originalPoint_2)));
		}

		TEST_METHOD(ArithmeticStructure_MatrixDegToRadFunctionTest)
		{
			float degToBeConvertedToRad{ 180 };
			float expectedRad{ M_PI };

			Assert::AreEqual(ArithmeticStructures::getRadiansForDeg(degToBeConvertedToRad), expectedRad);

			degToBeConvertedToRad*=  -1.0 ;
			expectedRad *= -1.0;

			Assert::AreEqual(ArithmeticStructures::getRadiansForDeg(degToBeConvertedToRad), expectedRad);

			degToBeConvertedToRad = 0.0;
			expectedRad = 0.0;
			Assert::AreEqual(ArithmeticStructures::getRadiansForDeg(degToBeConvertedToRad), expectedRad);
		}

		TEST_METHOD(ArithmeticStructure_MatrixRotationXAxisTest)
		{

			constexpr float rot_X_fullQuarter{M_PI_2}; // == 90°

			// is the transformation matrix created as expected
			const ArithmeticStructures::row4x4 m0_expected{ {1.0, 0.0, 0.0,0.0} }, m1_expected{ {0.0, cos(rot_X_fullQuarter), -1.0*sin(rot_X_fullQuarter), 0.0} },
				m2_expected{ {0.0,sin(rot_X_fullQuarter),cos(rot_X_fullQuarter), 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			const ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };
			auto transformationMatrix{ ArithmeticStructures::getRotationMatrix_XAxis(rot_X_fullQuarter) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does rotation work as expected - M * p = p`
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 0.0,1.0,0.0,1.0 };
			
			ArithmeticStructures::HomogenousCoordinates expectedResult_rotatedPoint{ 0.0, 0.0,1.0, 1.0  };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));

			constexpr float rot_X_halfQuarter{ M_PI_4 }; // == 45°
			constexpr float halfSqrt2{ M_SQRT2 * 0.5 };
			ArithmeticStructures::HomogenousCoordinates expectedResult_rotatedPoint_2{ 0.0, halfSqrt2,halfSqrt2, 1.0 };
			auto transformationMatrix_2{ ArithmeticStructures::getRotationMatrix_XAxis(rot_X_halfQuarter) };
			const auto rotatedPoint{ ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix_2, originalPoint) };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedPoint_2, rotatedPoint));
			
			// and also test whether the inversion works - M^-1 * p` = p
			ArithmeticStructures::inverseMatrix(transformationMatrix_2);
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(originalPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix_2, rotatedPoint)));

		}

		TEST_METHOD(ArithmeticStructure_MatrixRotationYAxisTest)
		{

			constexpr float rot_Y_fullQuarter{ M_PI_2 }; // == 90°

			// is the transformation matrix created as expected
			const ArithmeticStructures::row4x4 m0_expected{ {cos(rot_Y_fullQuarter), 0.0, sin(rot_Y_fullQuarter),0.0} }, m1_expected{ {0.0, 1.0,0.0, 0.0} },
				m2_expected{ {-1.0*sin(rot_Y_fullQuarter),0.0, cos(rot_Y_fullQuarter), 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			const ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };
			auto transformationMatrix{ ArithmeticStructures::getRotationMatrix_YAxis(rot_Y_fullQuarter) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does rotation work as expected - M * p = p`
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 0.0,0.0,1.0,1.0 };

			ArithmeticStructures::HomogenousCoordinates expectedResult_rotatedPoint{ 1.0, 0.0,0.0, 1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));

			constexpr float rot_Y_halfQuarter{ M_PI_4 }; // == 45°
			constexpr float halfSqrt2{ M_SQRT2 * 0.5 };
			ArithmeticStructures::HomogenousCoordinates expectedResult_rotatedPoint_2{ halfSqrt2,0.0,halfSqrt2, 1.0 };
			auto transformationMatrix_2{ ArithmeticStructures::getRotationMatrix_YAxis(rot_Y_halfQuarter) };
			const auto rotatedPoint{ ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix_2, originalPoint) };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedPoint_2, rotatedPoint));

			//// and also test whether the inversion works - M^-1 * p` = p
			ArithmeticStructures::inverseMatrix(transformationMatrix_2);
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(originalPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix_2, rotatedPoint)));

		}

		TEST_METHOD(ArithmeticStructure_MatrixRotationZAxisTest)
		{

			constexpr float rot_Z_fullQuarter{ M_PI_2 }; // == 90°

			// is the transformation matrix created as expected
			const ArithmeticStructures::row4x4 m0_expected{ {cos(rot_Z_fullQuarter), -1.0* sin(rot_Z_fullQuarter),0.0,0.0} }, m1_expected{ {sin(rot_Z_fullQuarter), cos(rot_Z_fullQuarter),0.0, 0.0} },
				m2_expected{ {0.0,0.0, 1.0, 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			const ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };
			auto transformationMatrix{ ArithmeticStructures::getRotationMatrix_ZAxis(rot_Z_fullQuarter) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does rotation work as expected - M * p = p`
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 0.0,1.0,0.0,1.0 };

			ArithmeticStructures::HomogenousCoordinates expectedResult_rotatedPoint{ -1.0, 0.0,0.0, 1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));

			constexpr float rot_Z_halfQuarter{ M_PI_4 }; // == 45°
			constexpr float halfSqrt2{ M_SQRT2 * 0.5 };
			ArithmeticStructures::HomogenousCoordinates expectedResult_rotatedPoint_2{ -1.0*halfSqrt2,halfSqrt2,0.0, 1.0 };
			auto transformationMatrix_2{ ArithmeticStructures::getRotationMatrix_ZAxis(rot_Z_halfQuarter) };
			const auto rotatedPoint{ ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix_2, originalPoint) };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedPoint_2, rotatedPoint));

			////// and also test whether the inversion works - M^-1 * p` = p
			ArithmeticStructures::inverseMatrix(transformationMatrix_2);
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(originalPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix_2, rotatedPoint)));

		}

		// xPropToY
		TEST_METHOD(ArithmeticStructure_MatrixShearingTest_1)
		{
			// shear (aka "move") x in proportion to y
			constexpr float propXtoY{ 1.0 };
			constexpr float propXtoZ{ 0.0 };
			constexpr float propYtoX{ 0.0 };
			constexpr float propYtoZ{ 0.0 };
			constexpr float propZtoX{ 0.0 };
			constexpr float propZtoY{ 0.0 };
			// is the transformation matrix created as expected
			
			const ArithmeticStructures::row4x4 m0_expected{ {1.0, propXtoY,propXtoZ,0.0} }, m1_expected{ {propYtoX, 1.0,propYtoZ, 0.0} },
				m2_expected{ {propZtoX,propZtoY, 1.0, 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			const ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };
			
			auto transformationMatrix{ ArithmeticStructures::getShearingMatrix(propXtoY, propXtoZ, propYtoX, propYtoZ, propZtoX, propZtoY) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does shearing work as expected - M * p = p`
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 2.0,3.0,4.0,1.0 };

			ArithmeticStructures::HomogenousCoordinates expectedResult_shearedPoint{ 5.0, 3.0,4.0, 1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_shearedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));
		}

		// xPropToZ
		TEST_METHOD(ArithmeticStructure_MatrixShearingTest_2)
		{
			// shear (aka "move") x in proportion to y
			constexpr float propXtoY{ 0.0 };
			constexpr float propXtoZ{ 1.0 };
			constexpr float propYtoX{ 0.0 };
			constexpr float propYtoZ{ 0.0 };
			constexpr float propZtoX{ 0.0 };
			constexpr float propZtoY{ 0.0 };
			// is the transformation matrix created as expected

			const ArithmeticStructures::row4x4 m0_expected{ {1.0, propXtoY,propXtoZ,0.0} }, m1_expected{ {propYtoX, 1.0,propYtoZ, 0.0} },
				m2_expected{ {propZtoX,propZtoY, 1.0, 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			const ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };

			auto transformationMatrix{ ArithmeticStructures::getShearingMatrix(propXtoY, propXtoZ, propYtoX, propYtoZ, propZtoX, propZtoY) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does shearing work as expected - M * p = p`
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 2.0,3.0,4.0,1.0 };

			ArithmeticStructures::HomogenousCoordinates expectedResult_shearedPoint{ 6.0, 3.0,4.0, 1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_shearedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));
		}

		//// yPropToX
		TEST_METHOD(ArithmeticStructure_MatrixShearingTest_3)
		{
			// shear (aka "move") x in proportion to y
			constexpr float propXtoY{ 0.0 };
			constexpr float propXtoZ{ 0.0 };
			constexpr float propYtoX{ 1.0 };
			constexpr float propYtoZ{ 0.0 };
			constexpr float propZtoX{ 0.0 };
			constexpr float propZtoY{ 0.0 };
			// is the transformation matrix created as expected

			const ArithmeticStructures::row4x4 m0_expected{ {1.0, propXtoY,propXtoZ,0.0} }, m1_expected{ {propYtoX, 1.0,propYtoZ, 0.0} },
				m2_expected{ {propZtoX,propZtoY, 1.0, 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			const ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };

			auto transformationMatrix{ ArithmeticStructures::getShearingMatrix(propXtoY, propXtoZ, propYtoX, propYtoZ, propZtoX, propZtoY) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does shearing work as expected - M * p = p`
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 2.0,3.0,4.0,1.0 };

			ArithmeticStructures::HomogenousCoordinates expectedResult_shearedPoint{ 2.0, 5.0,4.0, 1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_shearedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));
		}

		//// yPropToZ
		TEST_METHOD(ArithmeticStructure_MatrixShearingTest_4)
		{
			// shear (aka "move") x in proportion to y
			constexpr float propXtoY{ 0.0 };
			constexpr float propXtoZ{ 0.0 };
			constexpr float propYtoX{ 0.0 };
			constexpr float propYtoZ{ 1.0 };
			constexpr float propZtoX{ 0.0 };
			constexpr float propZtoY{ 0.0 };
			// is the transformation matrix created as expected

			const ArithmeticStructures::row4x4 m0_expected{ {1.0, propXtoY,propXtoZ,0.0} }, m1_expected{ {propYtoX, 1.0,propYtoZ, 0.0} },
				m2_expected{ {propZtoX,propZtoY, 1.0, 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			const ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };

			auto transformationMatrix{ ArithmeticStructures::getShearingMatrix(propXtoY, propXtoZ, propYtoX, propYtoZ, propZtoX, propZtoY) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does shearing work as expected - M * p = p`
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 2.0,3.0,4.0,1.0 };

			ArithmeticStructures::HomogenousCoordinates expectedResult_shearedPoint{ 2.0, 7.0,4.0, 1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_shearedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));
		}

		//// zPropToX
		TEST_METHOD(ArithmeticStructure_MatrixShearingTest_5)
		{
			// shear (aka "move") x in proportion to y
			constexpr float propXtoY{ 0.0 };
			constexpr float propXtoZ{ 0.0 };
			constexpr float propYtoX{ 0.0 };
			constexpr float propYtoZ{ 0.0 };
			constexpr float propZtoX{ 1.0 };
			constexpr float propZtoY{ 0.0 };
			// is the transformation matrix created as expected

			const ArithmeticStructures::row4x4 m0_expected{ {1.0, propXtoY,propXtoZ,0.0} }, m1_expected{ {propYtoX, 1.0,propYtoZ, 0.0} },
				m2_expected{ {propZtoX,propZtoY, 1.0, 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			const ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };

			auto transformationMatrix{ ArithmeticStructures::getShearingMatrix(propXtoY, propXtoZ, propYtoX, propYtoZ, propZtoX, propZtoY) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does shearing work as expected - M * p = p`
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 2.0,3.0,4.0,1.0 };

			ArithmeticStructures::HomogenousCoordinates expectedResult_shearedPoint{ 2.0, 3.0,6.0, 1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_shearedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));
		}

		//// zPropToY
		TEST_METHOD(ArithmeticStructure_MatrixShearingTest_6)
		{
			// shear (aka "move") x in proportion to y
			constexpr float propXtoY{ 0.0 };
			constexpr float propXtoZ{ 0.0 };
			constexpr float propYtoX{ 0.0 };
			constexpr float propYtoZ{ 0.0 };
			constexpr float propZtoX{ 0.0 };
			constexpr float propZtoY{ 1.0 };
			// is the transformation matrix created as expected

			const ArithmeticStructures::row4x4 m0_expected{ {1.0, propXtoY,propXtoZ,0.0} }, m1_expected{ {propYtoX, 1.0,propYtoZ, 0.0} },
				m2_expected{ {propZtoX,propZtoY, 1.0, 0.0} }, m3_expected{ {0.0,0.0,0.0,1.0} };
			const ArithmeticStructures::Matrix4x4 m_expected{ m0_expected, m1_expected, m2_expected, m3_expected };

			auto transformationMatrix{ ArithmeticStructures::getShearingMatrix(propXtoY, propXtoZ, propYtoX, propYtoZ, propZtoX, propZtoY) };

			Assert::IsTrue(ArithmeticStructures::matricesAreEqual_4x4(m_expected, transformationMatrix));

			// does shearing work as expected - M * p = p`
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 2.0,3.0,4.0,1.0 };

			ArithmeticStructures::HomogenousCoordinates expectedResult_shearedPoint{ 2.0, 3.0,7.0, 1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_shearedPoint, ArithmeticStructures::multiplyMatrixWithTuple(transformationMatrix, originalPoint)));
		}

		TEST_METHOD(ArithmeticStructure_MatrixConcatenationTest)
		{
			constexpr float rot_X_fullQuarter{ M_PI_2 }; // == 90°

			const auto rotationMatrix{ ArithmeticStructures::getRotationMatrix_XAxis(rot_X_fullQuarter) };

			// first apply rotation
			constexpr ArithmeticStructures::HomogenousCoordinates originalPoint{ 1.0,0.0,1.0,1.0 };

			ArithmeticStructures::HomogenousCoordinates expectedResult_rotatedPoint{ 1.0, -1.0,0.0, 1.0 };
			const auto rotatedPoint{ ArithmeticStructures::multiplyMatrixWithTuple(rotationMatrix, originalPoint) };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedPoint, rotatedPoint));


			// then apply scaling
			constexpr float scalingFactorXYZ{ 5.0 };

			const  auto scalingMatrix{ ArithmeticStructures::getScalingMatrix(scalingFactorXYZ, scalingFactorXYZ, scalingFactorXYZ) };

			ArithmeticStructures::HomogenousCoordinates expectedResult_rotatedAndScaledPoint{ 5.0, -5.0,0.0, 1.0 };
			const auto rotatedAndScaledPoint{ ArithmeticStructures::multiplyMatrixWithTuple(scalingMatrix, rotatedPoint) };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedAndScaledPoint, rotatedAndScaledPoint));

			// then apply translation
			constexpr float translation_X{ 10.0 }, translation_Y{ 5.0 }, translation_Z{ 7.0 };

			const auto translationMatrix{ ArithmeticStructures::getTranslationMatrix(translation_X, translation_Y, translation_Z)};

			ArithmeticStructures::HomogenousCoordinates expectedResult_rotatedAndScaledAndTranslatedPoint{ 15.0, 0.0,7.0, 1.0 };
			const auto rotatedAndScaledAndTranslatedPoint{ ArithmeticStructures::multiplyMatrixWithTuple(translationMatrix, rotatedAndScaledPoint) };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedAndScaledAndTranslatedPoint, rotatedAndScaledAndTranslatedPoint));

			const auto concatenatedMatrix{ArithmeticStructures::multiplyMatrices(ArithmeticStructures::multiplyMatrices(translationMatrix,scalingMatrix), rotationMatrix)};

			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedResult_rotatedAndScaledAndTranslatedPoint, ArithmeticStructures::multiplyMatrixWithTuple(concatenatedMatrix, originalPoint)));
		}

		TEST_METHOD(Ray_InitializationTest)
		{
			const ArithmeticStructures::HomogenousCoordinates origin{ 1.0,2.0,3.0,1.0 };
			const ArithmeticStructures::HomogenousCoordinates direction{ 4.0,5.0,6.0,0.0 };
			Ray ray{origin, direction};
			const ArithmeticStructures::HomogenousCoordinates expected_origin{ 1.0,2.0,3.0,1.0 };
			const ArithmeticStructures::HomogenousCoordinates expected_direction{ 4.0,5.0,6.0,0.0 };

			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expected_origin, ray.getOrigin()));
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expected_direction, ray.getDirection()));
		}

		TEST_METHOD(Ray_PositionTest)
		{
			constexpr float t_1{ 0.0 }, t_2{ 1.0 }, t_3{ -1.0 }, t_4{ 2.5 };

			const ArithmeticStructures::HomogenousCoordinates origin{ 2.0,3.0,4.0,1.0 };
			const ArithmeticStructures::HomogenousCoordinates direction{ 1.0,0.0,0.0,0.0 };
			Ray ray{ origin, direction };
			
			const ArithmeticStructures::HomogenousCoordinates expectedPos_t1{ 2.0,3.0,4.0,1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedPos_t1, ray.getPosition(t_1)));
			const ArithmeticStructures::HomogenousCoordinates expectedPos_t2{ 3.0,3.0,4.0,1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedPos_t2, ray.getPosition(t_2)));
			const ArithmeticStructures::HomogenousCoordinates expectedPos_t3{ 1.0,3.0,4.0,1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedPos_t3, ray.getPosition(t_3)));
			const ArithmeticStructures::HomogenousCoordinates expectedPos_t4{ 4.5,3.0,4.0,1.0 };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedPos_t4, ray.getPosition(t_4)));
		}

		TEST_METHOD(Ray_TranslationTest)
		{
			const ArithmeticStructures::HomogenousCoordinates origin{ 1.0,2.0,3.0,1.0 };
			const ArithmeticStructures::HomogenousCoordinates direction{ 0.0,1.0,0.0,0.0 };
			Ray original_ray{ origin, direction };
			constexpr float shift_x{ 3.0 }, shift_y{ 4.0 }, shift_z{ 5.0 };

			auto shiftedRay{ original_ray.translate(shift_x, shift_y, shift_z) };
			const auto expectedOrigin{ ArithmeticStructures::HomogenousCoordinates{4.0,6.0,8.0,1.0} };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedOrigin, shiftedRay.getOrigin()));
			const auto expectedDirection{ original_ray.getDirection() };
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedDirection, shiftedRay.getDirection()));
		}

		TEST_METHOD(SceneObject_SphereIntersectionTest)
		{
			// place the ray "in front" of the origin of the sphere and propagate the ray into the direction of the sphere, along the zaxis
			const ArithmeticStructures::HomogenousCoordinates ray_Origin{ 0.0,0.0,-5.0,1.0 };
			const ArithmeticStructures::HomogenousCoordinates ray_Direction{ 0.0,0.0,1.0,0.0 };
			Ray ray{ ray_Origin, ray_Direction };
			const ArithmeticStructures::HomogenousCoordinates sphere_Origin{ 0.0,0.0,0.0,1.0 };
			constexpr int sphere_Radius{ 1 };
			GeometricStructures::Sphere sphere{sphere_Origin, sphere_Radius};
			SceneObject sO{sphere};
			// expect the ray to intersect with the sphere two times. once on sphere entry, and afterwards while exiting the sphere
			SceneObject::Intersections expectedIntersections{ 4.0,6.0 };

			SceneObject::Intersections actualIntersections{ sO.getSphereIntersections(ray) };
			// first check whether all expected intersections have been found
			Assert::IsTrue(expectedIntersections.size() == actualIntersections.size());
			// then check whether they are correct
			Assert::AreEqual(expectedIntersections.at(0), actualIntersections.at(0), 0.0001f);
			Assert::AreEqual(expectedIntersections.at(1), actualIntersections.at(1), 0.0001f);


			// place the ray at the height of the sphere radius and intersect it tangentially
			ray.setOrigin(0.0, 1.0, - 5.0);
			
			// albeit only intersecting the speher tangentially, 2 (equal!) points of intersection will be returned
			expectedIntersections.at(0) = 5.0;
			expectedIntersections.at(1) = 5.0;

			actualIntersections = sO.getSphereIntersections(ray) ;
			// first check whether all expected intersections have been found
			Assert::IsTrue(expectedIntersections.size() == actualIntersections.size());
			// then check whether they are correct
			Assert::AreEqual(expectedIntersections.at(0), actualIntersections.at(0), 0.0001f);
			Assert::AreEqual(expectedIntersections.at(1), actualIntersections.at(1), 0.0001f);

			// place the ray above the sphere -> dont intersect the sphere at all
			ray.setOrigin(0.0, 2.0, -5.0);
			actualIntersections = sO.getSphereIntersections(ray);
			expectedIntersections.at(0) = SceneObject::Invalid;
			expectedIntersections.at(1) = SceneObject::Invalid;
			// no intersections expected
			Assert::AreEqual(expectedIntersections.at(0), actualIntersections.at(0), 0.0001f);
			Assert::AreEqual(expectedIntersections.at(1), actualIntersections.at(1), 0.0001f);

			// place the ray in the center of the sphere and expect 2 intersections. noe "behind" and one "in front" of the rays origin
			ray.setOrigin(0.0, 0.0, 0.0);
			expectedIntersections.at(0) = -1.0;
			expectedIntersections.at(1) = 1.0;

			actualIntersections = sO.getSphereIntersections(ray);
			// first check whether all expected intersections have been found
			Assert::IsTrue(expectedIntersections.size() == actualIntersections.size());
			// then check whether they are correct
			Assert::AreEqual(expectedIntersections.at(0), actualIntersections.at(0), 0.0001f);
			Assert::AreEqual(expectedIntersections.at(1), actualIntersections.at(1), 0.0001f);

			// place the ray "in front" of the sphere and the ray direction pointing "away" from the sphere. expecting 2 intersections "behind" the rays origin
			ray.setOrigin(0.0, 0.0, 5.0);
			expectedIntersections.at(0) = -6.0;
			expectedIntersections.at(1) = -4.0;

			actualIntersections = sO.getSphereIntersections(ray);
			// first check whether all expected intersections have been found
			Assert::IsTrue(expectedIntersections.size() == actualIntersections.size());
			// then check whether they are correct
			Assert::AreEqual(expectedIntersections.at(0), actualIntersections.at(0), 0.0001f);
			Assert::AreEqual(expectedIntersections.at(1), actualIntersections.at(1), 0.0001f);
		}

		TEST_METHOD(SceneObject_SphereHitTest)
		{
			// place the ray "in front" of the origin of the sphere and propagate the ray into the direction of the sphere, along the zaxis
			const ArithmeticStructures::HomogenousCoordinates ray_Origin{ 0.0,0.0,-5.0,1.0 };
			const ArithmeticStructures::HomogenousCoordinates ray_Direction{ 0.0,0.0,1.0,0.0 };
			Ray ray{ ray_Origin, ray_Direction };
			const ArithmeticStructures::HomogenousCoordinates sphere_Origin{ 0.0,0.0,0.0,1.0 };
			constexpr int sphere_Radius{ 1 };
			GeometricStructures::Sphere sphere{ sphere_Origin, sphere_Radius };
			SceneObject sO{ sphere };
			// expect the ray to intersect with the sphere two times. once on sphere entry, and afterwards while exiting the sphere
			// the entrypoint is expected to be the hit
			float expectedHit{ 4.0 };

			float actualHit{ sO.getSphereHit(ray) };
			
			Assert::AreEqual(expectedHit, actualHit, 0.0001f);


			// place the ray in the center of the sphere and expect 2 intersections. one "behind" and one "in front" of the rays origin
			// the expected hitpoint is the intersection "in front" of the ray
			ray.setOrigin(0.0, 0.0, 0.0);
			expectedHit = 1.0;

			actualHit =sO.getSphereHit(ray) ;
			Assert::AreEqual(expectedHit, actualHit, 0.0001f);

			// place the ray "in front" of the sphere and the ray direction pointing "away" from the sphere. expecting 2 intersections "behind" the rays origin
			// as there will be no hitpoint, the expected hit value is set to the according "indicator value" (NAN) -> comparing the actual value with isnan as NANs cant be compared directly
			ray.setOrigin(0.0, 0.0, 5.0);

			actualHit = sO.getSphereHit(ray);
			Assert::IsTrue(std::isnan(actualHit));
		}

		TEST_METHOD(Canvas_DimTest)
		{
			constexpr int xDim{ 256 }, yDim{ 354 };
			Canvas c(xDim, yDim);
			Assert::AreEqual(xDim, c.getDimX());
			Assert::AreEqual(yDim, c.getDimY());
		}

		TEST_METHOD(Canvas_SetImageDataTest)
		{
			constexpr int xDim{ 256 }, yDim{ 354 };
			Canvas c(xDim, yDim);
			ArithmeticStructures::HomogenousCoordinates imageData{ 1.0,0.0,0.0,1.0 };
			Assert::IsTrue(c.setImageData(12, 3, imageData));
			Assert::IsFalse(c.setImageData(345, 233, imageData));
		}

		TEST_METHOD(Canvas_GetImageDataTest)
		{
			constexpr int xDim{ 256 }, yDim{ 354 };
			constexpr int xCoord{ 12 }, yCoord{ 34 };
			const ArithmeticStructures::HomogenousCoordinates expectedColor{ 1.0,0.0,0.0,1.0 };
			Canvas c(xDim, yDim);
			c.setImageData(xCoord, yCoord, expectedColor);
			Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(expectedColor, c.getImageData(xCoord, yCoord)));
		}

		TEST_METHOD(PPMWriter_CreatePPMFile)
		{
			
			const std::string expectedFilename{ "unittest_filename.ppm" };
			// make sure a file with the expectedFilename is not present before executing the test
			std::ifstream istrm(expectedFilename, std::ios::binary);
			Assert::IsTrue(!istrm.is_open(), L"there shouldnt be a file with the expected filename before executing the test");
			istrm.close();

			PPMWriter imageWriter{ 256, 256, expectedFilename };
			imageWriter.createPPM(Canvas{ 0,0 });
			istrm.open(expectedFilename, std::ios::binary);
			Assert::IsTrue(istrm.is_open(), L"no file has been generated by the function under test");
			// todo: if test fails, close file again
			
			istrm.close();
			// cleanup folder after test
			std::filesystem::remove(expectedFilename);

			
		}

		TEST_METHOD(PPMWriter_CheckPPMFileMetaData)
		{
			const std::string createdFileName{ "createdFile.ppm" };
			// make sure a file with the expectedFilename is not present before executing the test
			std::ifstream createdFS(createdFileName, std::ios::binary);
			Assert::IsFalse(createdFS.is_open(), L"there shouldnt be a file with the filename before executing the test");
			createdFS.close();

			const std::string expectedPPMFormat{ "P3" };
			constexpr int expectedXDim{ 256 }, expectedYDim{ 256 };
			constexpr int expectedColorDepth{ 255 };

			PPMWriter imageWriter{ expectedXDim, expectedYDim, createdFileName };
			imageWriter.createPPM(Canvas{ 0,0 });

			std::string pPMFormatFromFile{ "" };
			int xDimFromFile{ 0 }, yDimFromFile{ 0 }, colorDepthFromFile{ 0 };

			createdFS.open(createdFileName, std::ios::binary);
			createdFS >> pPMFormatFromFile;
			createdFS >> xDimFromFile >> yDimFromFile;
			createdFS >> colorDepthFromFile;

			Assert::IsTrue(expectedPPMFormat == pPMFormatFromFile, L"read ppm format doesnt match expected ppm format");
			Assert::IsTrue(expectedXDim == xDimFromFile, L"read XDim doesnt match expected xDim");
			Assert::IsTrue(expectedYDim == yDimFromFile, L"read yDim doesnt match expected yDim");
			Assert::IsTrue(expectedColorDepth == colorDepthFromFile, L"read colorDepth doesnt match expected colorDepth");
			
			createdFS.close();
			// cleanup folder after test
			std::filesystem::remove(createdFileName);
		}

		TEST_METHOD(PPMWriter_CheckPPMFileImageData)
		{
			const std::string createdFileName{ "createdFile.ppm" };
			// make sure a file with the expectedFilename is not present before executing the test
			std::ifstream createdFS(createdFileName, std::ios::binary);
			Assert::IsFalse(createdFS.is_open(), L"there shouldnt be a file with the filename before executing the test");
			createdFS.close();

			constexpr int xDim{ 4 }, yDim{ 2 };

			PPMWriter imageWriter{ xDim, yDim, createdFileName };
			Canvas referenceCanvas(xDim, yDim);


			for (auto x = 0; x < xDim; x++)
			{
				for (auto y = 0; y < yDim; y++)
				{
					referenceCanvas.setImageData(x, y, ArithmeticStructures::HomogenousCoordinates{ (int)(x*255.0/xDim),(int)0.0,(int)0.0,1.0 });
				}
			}

			imageWriter.createPPM(referenceCanvas);

			std::string pPMFormatFromFile{ "" };
			int xDimFromFile{ 0 }, yDimFromFile{ 0 }, colorDepthFromFile{ 0 };

			createdFS.open(createdFileName, std::ios::binary);
			createdFS >> pPMFormatFromFile;
			createdFS >> xDimFromFile >> yDimFromFile;
			createdFS >> colorDepthFromFile;

			int r, g, b;
			std::array<std::array < std::array<int, 3>, 4>, 2> buffer{};

			for (auto y = 0; y <yDim; y++)
			{
				for (auto x = 0; x < xDim; x++)
				{
					createdFS >> r;
					createdFS >> g;
					createdFS >> b;
					buffer.at(y).at(x) = std::array<int, 3>{r, g, b};
				}
			}


			for (auto y = 0; y < yDim; y++)
			{
				for (auto x = 0; x < xDim; x++)
				{
					auto [v1, v2, v3] = buffer.at(y).at(x);
					
					Assert::IsTrue(ArithmeticStructures::coordinatesAreEqual(referenceCanvas.getImageData(x, y), ArithmeticStructures::HomogenousCoordinates{ v1,v2,v3,1.0 }), L"the value read from the file doesnt match the expected value");
				}
			}
			
			

			createdFS.close();

			// cleanup folder after test
			std::filesystem::remove(createdFileName);
		}

		/*TEST_METHOD(PPMWriter_CheckPPMDimensionDontMatchCanvasDimension)
		{
			Assert::Fail(L"test implementation needed, make sure a too large canvas cant be passed to the ppmwriter");
		}*/



		/*TEST_METHOD(PPMWriter_CheckPPMFileNameTooLong)
		{
			
			Assert::Fail(L"this test needs proper implementation");
		}*/
	};

	
}
