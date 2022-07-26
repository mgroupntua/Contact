using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.FEM.Structural.Line
{
	public class ContactNodeToSurface3D : IStructuralElementType
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
		private readonly double penaltyFactor;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[] DisplacementVector { get; set; }
		private double ContactArea { get; }

		public ContactNodeToSurface3D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier,
			double contactArea)
		{
			this.penaltyFactor = penaltyFactorMultiplier * youngModulus;
			this.DisplacementVector = new double[15];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
		}
		public ContactNodeToSurface3D(IReadOnlyList<INode> nodes, double penaltyFactor, double contactArea)
		{
			this.penaltyFactor = penaltyFactor;
			this.DisplacementVector = new double[15];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
		}
		public ContactNodeToSurface3D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplier, 1d)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactNodeToSurface3D(IReadOnlyList<INode> nodes, double penaltyFactor, IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactor, 1.0)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public CellType CellType { get; } = CellType.Unknown;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		private double[] XUpdatedVector()
		{
			var xVectorUpdated = new double[15];
			for (var i = 0; i < 5; i++)
			{
				xVectorUpdated[3 * i] = Nodes[i].X + DisplacementVector[3 * i];
				xVectorUpdated[3 * i + 1] = Nodes[i].Y + DisplacementVector[3 * i + 1];
				xVectorUpdated[3 * i + 2] = Nodes[i].Z + DisplacementVector[3 * i + 2];
			}
			return xVectorUpdated;
		}

		private Dictionary<int, double[]> SurfaceVectors(double[,] da1, double[,] da2)
		{
			var xUpdated = XUpdatedVector();
			return new Dictionary<int, double[]>()
			{
				{ 1, Matrix.CreateFromArray(da1).Multiply(xUpdated).Scale(-1.0) },
				{ 2, Matrix.CreateFromArray(da2).Multiply(xUpdated).Scale(-1.0) }
			};
		}

		private double[,] MetricTensor(Dictionary<int, double[]> masterSurfaceVectors) => new double[,]
		{
			{ masterSurfaceVectors[1].DotProduct(masterSurfaceVectors[1]), masterSurfaceVectors[1].DotProduct(masterSurfaceVectors[2]) },
			{ masterSurfaceVectors[2].DotProduct(masterSurfaceVectors[1]), masterSurfaceVectors[2].DotProduct(masterSurfaceVectors[2]) }
		};

		private double MetricTensorDet(double[,] m)
		{
			var detm = m[0, 0] * m[1, 1] - m[1, 0] * m[0, 1];
			return detm;
		}

		private Matrix InverseMetricTensor(double[,] m)
		{
			var detm = MetricTensorDet(m);
			var mInv = Matrix.CreateFromArray(new double[,]{ { m[1, 1], -m[0, 1] }, {-m[1,0], m[0,0] } }).Scale(1.0 / detm);			
			return mInv;
		}

		private double[] NormalVector(double[,] metricTensor, Dictionary<int, double[]> masterSurfaceVectors)
		{
			var n = (masterSurfaceVectors[1].CrossProduct(masterSurfaceVectors[2])).Scale(1.0 / (Math.Sqrt(MetricTensorDet(metricTensor))));
			return n;
		}

		private double[] Calculate_f(Dictionary<int, double[]> masterSurfaceVectors, double[,] aMatrix, double[] xUpdated) => new double[]
		{
			Vector.CreateFromArray(masterSurfaceVectors[1]).DotProduct(Matrix.CreateFromArray(aMatrix) * Vector.CreateFromArray(xUpdated)),
			Vector.CreateFromArray(masterSurfaceVectors[2]).DotProduct(Matrix.CreateFromArray(aMatrix) * Vector.CreateFromArray(xUpdated))
		};

		private double Calculate_e(double[,] aMatrix, double[] xUpdated)
		{
			var e = Vector.CreateFromArray(new double[]
			{
				0.25*(xUpdated[0] - xUpdated[3] + xUpdated[6] - xUpdated[9]),
				0.25*(xUpdated[1] - xUpdated[4] + xUpdated[7] - xUpdated[10]),
				0.25*(xUpdated[2] - xUpdated[5] + xUpdated[8] - xUpdated[11])
			}).DotProduct(Matrix.CreateFromArray(aMatrix) * Vector.CreateFromArray(xUpdated));
			return e;
		}

		private Vector CalculateDeltaKsi(double detm, double[,] mTensor, double[] fVector, double e)
		{
			var scalar = 1.0 / (detm - Math.Pow(e, 2) + 2.0 * e * mTensor[0, 1]);
			var matrix = new double[,]
			{
				{mTensor[1,1], e-mTensor[0,1] },
				{e-mTensor[1,0], mTensor[0,0] }
			};
			var deltaKsi = (Matrix.CreateFromArray(matrix) * Vector.CreateFromArray(fVector)).Scale(scalar);
			return deltaKsi;
		}

		private Vector Project(double[] ksiVectorInitial)
		{
			var maxIterations = 1000;
			var tol = Math.Pow(10.0, -4.0);
			var norm = new double();
			var ksiVector = Vector.CreateFromArray(ksiVectorInitial);
			var xUpdated = XUpdatedVector();
			for (var i = 1; i <= maxIterations; i++)
			{
				var aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
				var masterSurfaceVectors = SurfaceVectors(aMatrices.Item2, aMatrices.Item3);
				var f = Calculate_f(masterSurfaceVectors, aMatrices.Item1, xUpdated);
				var e = Calculate_e(aMatrices.Item1, xUpdated);
				var m = MetricTensor(masterSurfaceVectors);
				var detm = MetricTensorDet(m);
				var deltaKsi = CalculateDeltaKsi(detm, m, f, e);
				ksiVector += deltaKsi;
				norm = (deltaKsi).Norm2();
				if (norm <= tol)
				{
					break;
				}
			}
			if (norm > tol)
			{
				throw new Exception("CPP not found in current iterations");
			}
			else
			{
				return ksiVector;
			}

		}

		private Tuple<double[,], double[,], double[,]> CalculatePositionMatrix(double ksi1, double ksi2)
		{
			var N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
			var N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
			var N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
			var N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);

			var dN11 = -1.0 / 4.0 * (1.0 - ksi2);
			var dN21 = 1.0 / 4.0 * (1.0 - ksi2);
			var dN31 = 1.0 / 4.0 * (1.0 + ksi2);
			var dN41 = -1.0 / 4.0 * (1.0 + ksi2);

			var dN12 = -1.0 / 4.0 * (1.0 - ksi1);
			var dN22 = -1.0 / 4.0 * (1.0 + ksi1);
			var dN32 = 1.0 / 4.0 * (1.0 + ksi1);
			var dN42 = 1.0 / 4.0 * (1.0 - ksi1);

			var aMatrix = new double[,]
				{
					{ -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0, 0.0, 0.0 },
					{ 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0, 0.0 },
					{ 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0 }
				};

			var da1Matrix = new double[,]
				{
					{ -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 }
				};

			var da2Matrix = new double[,]
				{
					{ -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 }
				};
			return new Tuple<double[,], double[,], double[,]>(aMatrix, da1Matrix, da2Matrix);
		}

		private double CalculatePenetration(double[] normalVector, double[,] aMatrix, double[] xUpdated)
		{
			var ksi3 = Vector.CreateFromArray(xUpdated).DotProduct(Matrix.CreateFromArray(aMatrix).Transpose() * Vector.CreateFromArray(normalVector));
			return ksi3;
		}
		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		private Matrix CalculateMainStiffnessPart(double[] n, Matrix aMatrix)
		{
			var nxn = n.TensorProduct(n);
			var nxn_A = nxn.MultiplyRight(aMatrix);
			var AT_nxn_A = aMatrix.Transpose().MultiplyRight(nxn_A);
			var mainStiffnessMatrix = AT_nxn_A.Scale(penaltyFactor);
			return mainStiffnessMatrix;
		}

		private Matrix CalculateRotationalStiffnessPart(double[] normalVector, Matrix aMatrix, Matrix dAd1Matrix, Matrix dAd2Matrix, Dictionary<int, double[]> dRho, double ksi3, double[,] metricTensor)
		{
			var mInv = InverseMetricTensor(metricTensor);
			var scalar1 = penaltyFactor * ksi3 * mInv[0, 0];
			var scalar2 = penaltyFactor * ksi3 * mInv[1, 0];
			var scalar3 = penaltyFactor * ksi3 * mInv[0, 1];
			var scalar4 = penaltyFactor * ksi3 * mInv[1, 1];
			var mat1 = dAd1Matrix.Transpose() * (normalVector.TensorProduct(dRho[1])) * aMatrix + aMatrix.Transpose() * (dRho[1].TensorProduct(normalVector)) * dAd1Matrix;
			var mat2 = dAd1Matrix.Transpose() * (normalVector.TensorProduct(dRho[2])) * aMatrix + aMatrix.Transpose() * (dRho[1].TensorProduct(normalVector)) * dAd2Matrix;
			var mat3 = dAd2Matrix.Transpose() * (normalVector.TensorProduct(dRho[1])) * aMatrix + aMatrix.Transpose() * (dRho[2].TensorProduct(normalVector)) * dAd1Matrix;
			var mat4 = dAd2Matrix.Transpose() * (normalVector.TensorProduct(dRho[2])) * aMatrix + aMatrix.Transpose() * (dRho[2].TensorProduct(normalVector)) * dAd2Matrix;
			var Kr = mat1.Scale(scalar1) + mat2.Scale(scalar2) + mat3.Scale(scalar3) + mat4.Scale(scalar4);
			return Kr;
		}

		public IMatrix StiffnessMatrix()
		{
			var ksiVector = Project(new double[2]);
			if (Math.Abs(ksiVector[0]) <= 1.05 && Math.Abs(ksiVector[1]) <= 1.05)
			{
				var aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
				var masterSurfaceVectors = SurfaceVectors(aMatrices.Item2, aMatrices.Item3);
				var m = MetricTensor(masterSurfaceVectors);
				var n = NormalVector(m, masterSurfaceVectors);
				var xUpdated = XUpdatedVector();
				var ksi3 = CalculatePenetration(n, aMatrices.Item1, xUpdated);
				if (ksi3 <= 0)
				{
					var Km = CalculateMainStiffnessPart(n,Matrix.CreateFromArray(aMatrices.Item1));
					var Kr = CalculateRotationalStiffnessPart(n, Matrix.CreateFromArray(aMatrices.Item1),
						Matrix.CreateFromArray(aMatrices.Item2), Matrix.CreateFromArray(aMatrices.Item3),
						masterSurfaceVectors, ksi3, m);
					var K = Km.Add(Kr);
					return dofEnumerator.GetTransformedMatrix(K);
				}
				else
				{
					return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(new double[15, 15]));
				}
			}
			else
			{
				return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(new double[15, 15]));
			}
		}
		public IMatrix PhysicsMatrix() => StiffnessMatrix();
		public IMatrix MassMatrix()
		{
			double[,] massMatrix = new double[15, 15];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(massMatrix));
		}

		public IMatrix DampingMatrix()
		{
			double[,] dampingMatrix = new double[15, 15];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(dampingMatrix));
		}

		private double[] CreateInternalGlobalForcesVector()
		{
			var ksiVector = Project(new double[2]);
			if (Math.Abs(ksiVector[0]) <= 1.05 && Math.Abs(ksiVector[1]) <= 1.05)
			{
				var aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
				var masterSurfaceVectors = SurfaceVectors(aMatrices.Item2, aMatrices.Item3);
				var m = MetricTensor(masterSurfaceVectors);
				var n = NormalVector(m, masterSurfaceVectors);
				var xUpdated = XUpdatedVector();
				var ksi3 = CalculatePenetration(n, aMatrices.Item1, xUpdated);
				if (ksi3 <= 0)
				{
					var AT = Matrix.CreateFromArray(aMatrices.Item1).Transpose();
					var internalGlobalForcesVector = AT.Multiply(n).Scale(penaltyFactor * ksi3);
					return internalGlobalForcesVector;
				}
				else
				{
					var internalGlobalForcesVector = new double[15];
					return internalGlobalForcesVector;
				}
			}
			else
			{
				var internalGlobalForcesVector = new double[15];
				return internalGlobalForcesVector;
			}
		}

		public Tuple<double[], double[]> CalculateResponse(double[] local_Displacements)
		{
			// WARNING: 1) No strains are computed 2) localdDisplacements are not used.
			double[] strains = null;
			double[] forces = CreateInternalGlobalForcesVector();
			double[] stresses = Array.ConvertAll(forces, x => x / ContactArea);
			if (DisplacementVector == null || DisplacementVector.Length != local_Displacements.Length)
			{
				DisplacementVector = new double[local_Displacements.Length];
			}

			Array.Copy(local_Displacements, DisplacementVector, local_Displacements.Length);

			return new Tuple<double[], double[]>(strains, stresses);
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			return CalculateResponseIntegral();
		}

		public double[] CalculateResponseIntegral() => CreateInternalGlobalForcesVector();

		public void SaveConstitutiveLawState() { }

		#endregion

		#region IFiniteElement Members


		public bool ConstitutiveLawModified => false;

		public void ResetConstitutiveLawModified() { }

		#endregion

		#region IFiniteElement Members

		public void ClearConstitutiveLawState() { }

		public void ClearConstitutiveLawStresses() => throw new NotImplementedException();

		#endregion
	}
}
