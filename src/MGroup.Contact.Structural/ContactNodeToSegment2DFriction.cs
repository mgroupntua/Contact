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
	public class ContactNodeToSegment2DFriction : IStructuralElementType
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
		private readonly double PenaltyFactorNormal;
		private readonly double PenaltyFactorTangential;
		private readonly double StickingCoefficient;
		private readonly double SlidingCoefficient;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[] DisplacementVector { get; set; }
		private double ContactArea { get; }

		public ContactNodeToSegment2DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient,
			double contactArea)
		{
			this.PenaltyFactorNormal = penaltyFactorMultiplierNormal * youngModulus;
			this.PenaltyFactorTangential = penaltyFactorMultiplierTangential * youngModulus;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[6];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
		}
		public ContactNodeToSegment2DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal, double penaltyFactorTangential,
			double stickingCoefficient, double slidingCoefficient, double contactArea)
		{
			this.PenaltyFactorNormal = penaltyFactorNormal;
			this.PenaltyFactorTangential = penaltyFactorTangential;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[6];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
		}
		public ContactNodeToSegment2DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient,
			double contactArea, IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplierNormal, penaltyFactorMultiplierTangential,
					stickingCoefficient, slidingCoefficient, 1d)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactNodeToSegment2DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal, double penaltyFactorTangential,
			double stickingCoefficient, double slidingCoefficient, double contactArea, IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactorNormal, penaltyFactorTangential, stickingCoefficient, slidingCoefficient, 1.0)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public CellType CellType { get; } = CellType.Unknown;// I guess?

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		public double InitialProjectionPoint()
		{
			var Xm1Initial = Nodes[0].X;
			var Ym1Initial = Nodes[0].Y;
			var Xm2Initial = Nodes[1].X;
			var Ym2Initial = Nodes[1].Y;
			var XsInitial = Nodes[2].X;
			var YsInitial = Nodes[2].Y;
			var ksi1Initial = (2.0 * (XsInitial * (Xm2Initial - Xm1Initial) + YsInitial * (Ym2Initial - Ym1Initial)) - Math.Pow(Xm2Initial, 2) - Math.Pow(Ym2Initial, 2) + Math.Pow(Xm1Initial, 2) + Math.Pow(Ym1Initial, 2)) / (Math.Pow(Xm2Initial - Xm1Initial, 2) + Math.Pow(Ym2Initial - Ym1Initial, 2));
			return ksi1Initial;
		}

		public double ClosestPointProjection()
		{
			var Xm1 = Nodes[0].X + DisplacementVector[0];
			var Ym1 = Nodes[0].Y + DisplacementVector[1];
			var Xm2 = Nodes[1].X + DisplacementVector[2];
			var Ym2 = Nodes[1].Y + DisplacementVector[3];
			var Xs = Nodes[2].X + DisplacementVector[4];
			var Ys = Nodes[2].Y + DisplacementVector[5];

			var ksi1 = (2.0 * (Xs * (Xm2 - Xm1) + Ys * (Ym2 - Ym1)) - Math.Pow(Xm2, 2) - Math.Pow(Ym2, 2) + Math.Pow(Xm1, 2) + Math.Pow(Ym1, 2)) / (Math.Pow(Xm2 - Xm1, 2) + Math.Pow(Ym2 - Ym1, 2));
			return ksi1;
		}

		private double CalculateNormalGap(double[,] aMatrix, double[] n)
		{
			var AT = Matrix.CreateFromArray(aMatrix).Transpose();
			var AT_n = Vector.CreateFromArray(AT.Multiply(n));
			var xupd = Vector.CreateFromArray(new double[]
			{
				Nodes[0].X + DisplacementVector[0],
				Nodes[0].Y + DisplacementVector[1],
				Nodes[1].X + DisplacementVector[2],
				Nodes[1].Y + DisplacementVector[3],
				Nodes[2].X + DisplacementVector[4],
				Nodes[2].Y + DisplacementVector[5]
			});
			var normalGap = xupd.DotProduct(AT_n);
			return normalGap;
		}

		private Tuple<double[], double, double[], double, double[]> MasterSegmentGeoParameters(double[,] daMatrix)
		{
			var Xm1 = Nodes[0].X + DisplacementVector[0];
			var Ym1 = Nodes[0].Y + DisplacementVector[1];
			var Xm2 = Nodes[1].X + DisplacementVector[2];
			var Ym2 = Nodes[1].Y + DisplacementVector[3];
			var Xs = Nodes[2].X + DisplacementVector[4];
			var Ys = Nodes[2].Y + DisplacementVector[5];
			var xupd = new double[] { -Xm1, -Ym1, -Xm2, -Ym2, -Xs, -Ys };
			var surfaceVector = Matrix.CreateFromArray(daMatrix).Multiply(xupd);
			var detm = surfaceVector.DotProduct(surfaceVector);
			var m11 = 1.0 / detm;
			var vector = new double[] { Ym2 - Ym1, Xm1 - Xm2 };
			var scalarCoef = -1.0 / (2.0 * Math.Sqrt(detm));
			var normalUnitVec = vector.Scale(scalarCoef);
			var tVector = surfaceVector.Scale(1.0 / Math.Sqrt(detm));
			return new Tuple<double[], double, double[], double, double[]>(surfaceVector, m11, normalUnitVec, detm, tVector);
		}

		private Tuple<double[,], double[,]> CalculatePositionMatrix(double ksi1)
		{
			var N1 = 1.0 / 2.0 * (1.0 - ksi1);
			var N2 = 1.0 / 2.0 * (1.0 + ksi1);
			var dN1 = -1.0 / 2.0;
			var dN2 = 1.0 / 2.0;
			var aMatrix = new double[,]
				{
					{ -N1 ,0.0 ,-N2 ,0.0 ,1.0 ,0.0 },
					{0.0, -N1 , 0.0 ,-N2, 0.0, 1.0 }
				};

			var daMatrix = new double[,]
				{
					{ -dN1 ,0.0 ,-dN2 ,0.0 ,0.0 ,0.0 },
					{0.0, -dN1 , 0.0 ,-dN2, 0.0, 0.0 }
				};
			return new Tuple<double[,], double[,]>(aMatrix, daMatrix);
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		private Matrix CalculateMainStiffnessPart(double[] n, Matrix aMatrix)
		{
			var nxn = n.TensorProduct(n);
			var nxn_A = nxn.MultiplyRight(aMatrix);
			var AT_nxn_A = aMatrix.Transpose().MultiplyRight(nxn_A);
			var mainStiffnessMatrix = AT_nxn_A.Scale(PenaltyFactorNormal);
			return mainStiffnessMatrix;
		}

		private Matrix CalculateRotationalStiffnessPart(Matrix A, Matrix dA, double[] n, double penetration, double m11, double[] surfaceVector)
		{
			var coef = PenaltyFactorNormal * penetration * m11;
			var n_x_dRho = n.TensorProduct(surfaceVector);
			var dRho_x_n = surfaceVector.TensorProduct(n);
			var firstTerm = dA.Transpose().MultiplyRight(n_x_dRho.MultiplyRight(A));
			var secondTerm = A.Transpose().MultiplyRight(dRho_x_n.MultiplyRight(dA));
			var rotationalPart = (firstTerm + secondTerm).Scale(coef);
			return rotationalPart;
		}

		public IMatrix StiffnessMatrix()
		{
			var Ksi1Current = ClosestPointProjection();
			if (Math.Abs(Ksi1Current) <= 1.05)
			{
				var positionMatrices = CalculatePositionMatrix(Ksi1Current);
				var A = positionMatrices.Item1;
				var dA = positionMatrices.Item2;
				var masterSegmentGeometry = MasterSegmentGeoParameters(dA);
				var surfaceVector = masterSegmentGeometry.Item1;
				var metricTensorContravariant = masterSegmentGeometry.Item2;
				var n = masterSegmentGeometry.Item3;
				var detm = masterSegmentGeometry.Item4;
				var tVector = masterSegmentGeometry.Item5;
				var penetration = CalculateNormalGap(A, n);
				if (penetration <= 0)
				{
					var mainPart = CalculateMainStiffnessPart(n, Matrix.CreateFromArray(A));
					var rotationalPart = CalculateRotationalStiffnessPart(Matrix.CreateFromArray(A), Matrix.CreateFromArray(dA),
										n, penetration, metricTensorContravariant, surfaceVector);
					var globalStiffnessMatrixNormalPart = mainPart.Add(rotationalPart);
					var Ksi1Initial = InitialProjectionPoint();
					var deltaKsi = Ksi1Current - Ksi1Initial;
					var trialTangentialTraction = -PenaltyFactorTangential * detm * deltaKsi;
					var phi = Math.Sqrt(trialTangentialTraction * trialTangentialTraction * metricTensorContravariant) -
						StickingCoefficient * PenaltyFactorNormal * Math.Abs(penetration);
					var globalStiffnessMatrixTangentialPart = Matrix.CreateFromArray(new double[6, 6]);
					if (phi <= 0.0)
					{
						//Sticking
						var stickingK1 = (Matrix.CreateFromArray(A).Transpose() * tVector.TensorProduct(tVector) * Matrix.CreateFromArray(A)).Scale(PenaltyFactorTangential);
						var StickingK2 = (Matrix.CreateFromArray(dA).Transpose() * tVector.TensorProduct(tVector) * Matrix.CreateFromArray(A)).Scale(trialTangentialTraction * metricTensorContravariant);
						globalStiffnessMatrixTangentialPart = StickingK2 + StickingK2.Transpose() - stickingK1;
					}
					else
					{
						//Sliding
						var slidingK1 = (Matrix.CreateFromArray(A).Transpose() * tVector.TensorProduct(n) * Matrix.CreateFromArray(A)).Scale(SlidingCoefficient * PenaltyFactorNormal * (trialTangentialTraction / Math.Abs(trialTangentialTraction)));
						var slidingK2 = (Matrix.CreateFromArray(dA).Transpose() * tVector.TensorProduct(tVector) * Matrix.CreateFromArray(A)).
							Scale(SlidingCoefficient * PenaltyFactorNormal * Math.Abs(penetration) * Math.Sqrt(metricTensorContravariant) * 
							(trialTangentialTraction / Math.Abs(trialTangentialTraction)));
						globalStiffnessMatrixTangentialPart = slidingK2 + slidingK2.Transpose() - slidingK1;
					}
					var globalStifnessMatrix = globalStiffnessMatrixNormalPart.Add(globalStiffnessMatrixTangentialPart);
					return dofEnumerator.GetTransformedMatrix(globalStifnessMatrix);
				}
				else
				{
					var globalStifnessMatrix = new double[6, 6];
					return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(globalStifnessMatrix));
				}

			}
			else
			{
				var globalStifnessMatrix = new double[6, 6];
				return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(globalStifnessMatrix));
			}
		}

		public IMatrix PhysicsMatrix() => StiffnessMatrix();

		public IMatrix MassMatrix()
		{
			double[,] massMatrix = new double[6, 6];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(massMatrix));
		}

		public IMatrix DampingMatrix()
		{
			double[,] dampingMatrix = new double[6, 6];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(dampingMatrix));
		}

		private double[] CreateInternalGlobalForcesVector()
		{
			var Ksi1Current = ClosestPointProjection();
			if (Math.Abs(Ksi1Current) <= 1.05)
			{
				var positionMatrices = CalculatePositionMatrix(Ksi1Current);
				var A = positionMatrices.Item1;
				var dA = positionMatrices.Item2;
				var masterSegmentGeometry = MasterSegmentGeoParameters(dA);
				var metricTensorContravariant = masterSegmentGeometry.Item2;
				var nVector = masterSegmentGeometry.Item3;
				var detm = masterSegmentGeometry.Item4;
				var tVector = masterSegmentGeometry.Item5;
				var penetration = CalculateNormalGap(A, nVector);
				if (penetration <= 0)
				{
					var AT = Matrix.CreateFromArray(A).Transpose();
					var AT_n = AT.Multiply(nVector);
					var AT_t = AT.Multiply(tVector);
					var Ksi1Initial = InitialProjectionPoint();
					var deltaKsi = Ksi1Current - Ksi1Initial;
					var trialTangentialTraction = -PenaltyFactorTangential * detm * deltaKsi;
					var phi = Math.Sqrt(trialTangentialTraction * trialTangentialTraction * metricTensorContravariant) -
						StickingCoefficient * PenaltyFactorNormal * Math.Abs(penetration);
					if (phi <= 0.0)
					{
						var internalGlobalForcesVector = AT_n.Scale(PenaltyFactorNormal * penetration).Add
							(AT_t.Scale(trialTangentialTraction * Math.Sqrt(metricTensorContravariant)));
						return internalGlobalForcesVector;
					}
					else
					{
						var T1 = (trialTangentialTraction / Math.Abs(trialTangentialTraction)) * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(penetration) * Math.Sqrt(detm);
						var internalGlobalForcesVector = AT_n.Scale(PenaltyFactorNormal * penetration).Add
							(AT_t.Scale(T1 * Math.Sqrt(metricTensorContravariant)));
						return internalGlobalForcesVector;
					}
				}
				else
				{
					var internalGlobalForcesVector = new double[6];
					return internalGlobalForcesVector;
				}
			}
			else
			{
				var internalGlobalForcesVector = new double[6];
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

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements) => CalculateResponseIntegral();

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
