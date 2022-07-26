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
	public class ContactSegmentToSegment2DFriction : IStructuralElementType
	{
		private readonly IDofType[][] dofs;
		private readonly double PenaltyFactorNormal;
		private readonly double PenaltyFactorTangential;
		private readonly double StickingCoefficient;
		private readonly double SlidingCoefficient;
		private readonly int MasterSegmentOrder;
		private readonly int SlaveSegmentOrder;
		private readonly int IntegrationPoints;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[] DisplacementVector { get; set; }
		private double ContactArea { get; }
		private Dictionary<int, double> IntegrationPointsStickingPoints { get; set; }
		private Dictionary<int, double> IntegrationPointsTangentialTractions { get; set; }
		public ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient,
			double contactArea)
		{
			if (nodes.Count != 4)
			{
				throw new ArgumentException("This Constructor can be used only for linear segments");
			}
			this.PenaltyFactorNormal = penaltyFactorMultiplierNormal * youngModulus;
			this.PenaltyFactorTangential = penaltyFactorMultiplierTangential * youngModulus;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[2 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
						StructuralDof.TranslationX, StructuralDof.TranslationY
				};
			}
			this.MasterSegmentOrder = 1;
			this.SlaveSegmentOrder = 1;
			this.IntegrationPoints = 2;
			this.IntegrationPointsStickingPoints = new Dictionary<int, double>() 
			{
				{ 1, 0.0 },
				{ 2, 0.0 }
			};
			this.IntegrationPointsTangentialTractions = new Dictionary<int, double>()
			{
				{ 1, 0.0 },
				{ 2, 0.0 }
			};
		}
		public ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal, double penaltyFactorTangential,
			double stickingCoefficient, double slidingCoefficient, double contactArea)
		{
			if (nodes.Count != 4)
			{
				throw new ArgumentException("This Constructor can be used only for linear segments");
			}
			this.PenaltyFactorNormal = penaltyFactorNormal;
			this.PenaltyFactorTangential = penaltyFactorTangential;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[2 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
						StructuralDof.TranslationX, StructuralDof.TranslationY
				};
			}
			this.MasterSegmentOrder = 1;
			this.SlaveSegmentOrder = 1;
			this.IntegrationPoints = 2;
			this.IntegrationPointsStickingPoints = new Dictionary<int, double>()
			{
				{ 1, 0.0 },
				{ 2, 0.0 }
			};
			this.IntegrationPointsTangentialTractions = new Dictionary<int, double>()
			{
				{ 1, 0.0 },
				{ 2, 0.0 }
			};
		}
		public ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient, double contactArea,
			int masterSegmentOrder, int slaveSegmentOrder)
		{
			if ((masterSegmentOrder != 1 && masterSegmentOrder != 2) || (slaveSegmentOrder != 1 && slaveSegmentOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic segments can be defined");
			}
			if (nodes.Count != slaveSegmentOrder + masterSegmentOrder + 2)
			{
				throw new ArgumentException("Inconsistent input regarding the nodes & the order of the segments");
			}
			this.PenaltyFactorNormal = penaltyFactorMultiplierNormal * youngModulus;
			this.PenaltyFactorTangential = penaltyFactorMultiplierTangential * youngModulus;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[2 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
						StructuralDof.TranslationX, StructuralDof.TranslationY
				};
			}
			this.MasterSegmentOrder = masterSegmentOrder;
			this.SlaveSegmentOrder = slaveSegmentOrder;
			if (slaveSegmentOrder == 1)
			{
				this.IntegrationPoints = 2;
			}
			else
			{
				this.IntegrationPoints = 3;
			}
			var integrationPointsStickingPoints = new Dictionary<int, double>();
			var integrationPointsTangentialTractions = new Dictionary<int, double>();
			for (var i = 0; i < this.IntegrationPoints; i++)
			{
				integrationPointsStickingPoints.Add(i + 1, 0.0);
				integrationPointsTangentialTractions.Add(i + 1, 0.0);

			}
			this.IntegrationPointsStickingPoints = integrationPointsStickingPoints;
			this.IntegrationPointsTangentialTractions = integrationPointsTangentialTractions;
		}
		public ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal, double penaltyFactorTangential,
			double stickingCoefficient, double slidingCoefficient, double contactArea,
			int masterSegmentOrder, int slaveSegmentOrder)
		{
			if ((masterSegmentOrder != 1 && masterSegmentOrder != 2) || (slaveSegmentOrder != 1 && slaveSegmentOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic segments can be defined");
			}
			if (nodes.Count != slaveSegmentOrder + masterSegmentOrder + 2)
			{
				throw new ArgumentException("Inconsistent input regarding the nodes & the order of the segments");
			}
			this.PenaltyFactorNormal = penaltyFactorNormal;
			this.PenaltyFactorTangential = penaltyFactorTangential;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[2 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
						StructuralDof.TranslationX, StructuralDof.TranslationY
				};
			}
			this.MasterSegmentOrder = masterSegmentOrder;
			this.SlaveSegmentOrder = slaveSegmentOrder;
			if (slaveSegmentOrder == 1)
			{
				this.IntegrationPoints = 2;
			}
			else
			{
				this.IntegrationPoints = 3;
			}
			var integrationPointsStickingPoints = new Dictionary<int, double>();
			var integrationPointsTangentialTractions = new Dictionary<int, double>();
			for (var i = 0; i < this.IntegrationPoints; i++)
			{
				integrationPointsStickingPoints.Add(i + 1, 0.0);
				integrationPointsTangentialTractions.Add(i + 1, 0.0);

			}
			this.IntegrationPointsStickingPoints = integrationPointsStickingPoints;
			this.IntegrationPointsTangentialTractions = integrationPointsTangentialTractions;
		}
		public ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient, double contactArea,
			int masterSegmentOrder, int slaveSegmentOrder, int integrationPoints)
		{
			if ((masterSegmentOrder != 1 && masterSegmentOrder != 2) || (slaveSegmentOrder != 1 && slaveSegmentOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic segments can be defined");
			}
			if (nodes.Count != slaveSegmentOrder + masterSegmentOrder + 2)
			{
				throw new ArgumentException("Inconsistent input regarding the nodes & the order of the segments");
			}
			if (integrationPoints < 2 || integrationPoints > 10)
			{
				throw new ArgumentException("Between [2,10] Gauss points can be defined");
			}
			else if(slaveSegmentOrder == 2 && integrationPoints < 3)
			{
				throw new ArgumentException("insufficient integration order");
			}
			this.PenaltyFactorNormal = penaltyFactorMultiplierNormal * youngModulus;
			this.PenaltyFactorTangential = penaltyFactorMultiplierTangential * youngModulus;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[2 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
						StructuralDof.TranslationX, StructuralDof.TranslationY
				};
			}
			this.MasterSegmentOrder = masterSegmentOrder;
			this.SlaveSegmentOrder = slaveSegmentOrder;
			this.IntegrationPoints = integrationPoints;
			var integrationPointsStickingPoints = new Dictionary<int, double>();
			var integrationPointsTangentialTractions = new Dictionary<int, double>();
			for (var i = 0; i < this.IntegrationPoints; i++)
			{
				integrationPointsStickingPoints.Add(i + 1, 0.0);
				integrationPointsTangentialTractions.Add(i + 1, 0.0);

			}
			this.IntegrationPointsStickingPoints = integrationPointsStickingPoints;
			this.IntegrationPointsTangentialTractions = integrationPointsTangentialTractions;
		}
		public ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal, double penaltyFactorTangential,
			double stickingCoefficient, double slidingCoefficient, double contactArea,
			int masterSegmentOrder, int slaveSegmentOrder, int integrationPoints)
		{
			if ((masterSegmentOrder != 1 && masterSegmentOrder != 2) || (slaveSegmentOrder != 1 && slaveSegmentOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic segments can be defined");
			}
			if (nodes.Count != slaveSegmentOrder + masterSegmentOrder + 2)
			{
				throw new ArgumentException("Inconsistent input regarding the nodes & the order of the segments");
			}
			if (integrationPoints < 2 || integrationPoints > 10)
			{
				throw new ArgumentException("Between [2,10] Gauss points can be defined");
			}
			else if (slaveSegmentOrder == 2 && integrationPoints < 3)
			{
				throw new ArgumentException("insufficient integration order");
			}
			this.PenaltyFactorNormal = penaltyFactorNormal;
			this.PenaltyFactorTangential = penaltyFactorTangential;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[2 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
						StructuralDof.TranslationX, StructuralDof.TranslationY
				};
			}
			this.MasterSegmentOrder = masterSegmentOrder;
			this.SlaveSegmentOrder = slaveSegmentOrder;
			this.IntegrationPoints = integrationPoints;
			var integrationPointsStickingPoints = new Dictionary<int, double>();
			var integrationPointsTangentialTractions = new Dictionary<int, double>();
			for (var i = 0; i < this.IntegrationPoints; i++)
			{
				integrationPointsStickingPoints.Add(i + 1, 0.0);
				integrationPointsTangentialTractions.Add(i + 1, 0.0);

			}
			this.IntegrationPointsStickingPoints = integrationPointsStickingPoints;
			this.IntegrationPointsTangentialTractions = integrationPointsTangentialTractions;
		}
		ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient, double contactArea,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplierNormal, penaltyFactorMultiplierTangential,
					stickingCoefficient, slidingCoefficient, contactArea)
		{
			this.dofEnumerator = dofEnumerator;
		}
		ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal, double penaltyFactorTangential,
			double stickingCoefficient, double slidingCoefficient, double contactArea,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactorNormal, penaltyFactorTangential, stickingCoefficient,
				  slidingCoefficient, contactArea)
		{
			this.dofEnumerator = dofEnumerator;
		}
		ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient, double contactArea,
			int masterSegmentOrder, int slaveSegmentOrder,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplierNormal, penaltyFactorMultiplierTangential,
			stickingCoefficient, slidingCoefficient, contactArea, masterSegmentOrder, slaveSegmentOrder)
		{
			this.dofEnumerator = dofEnumerator;
		}
		ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal, double penaltyFactorTangential,
			double stickingCoefficient, double slidingCoefficient, double contactArea,
			int masterSegmentOrder, int slaveSegmentOrder,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactorNormal, penaltyFactorTangential, stickingCoefficient,
				  slidingCoefficient, contactArea, masterSegmentOrder, slaveSegmentOrder)
		{
			this.dofEnumerator = dofEnumerator;
		}
		ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient, double contactArea,
			int masterSegmentOrder, int slaveSegmentOrder, int integrationPoints,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplierNormal, penaltyFactorMultiplierTangential,
					stickingCoefficient, slidingCoefficient, contactArea, masterSegmentOrder, slaveSegmentOrder, integrationPoints)
		{
			this.dofEnumerator = dofEnumerator;
		}
		ContactSegmentToSegment2DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal, double penaltyFactorTangential,
			double stickingCoefficient, double slidingCoefficient, double contactArea,
			int masterSegmentOrder, int slaveSegmentOrder, int integrationPoints,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactorNormal, penaltyFactorTangential, stickingCoefficient,
				  slidingCoefficient, contactArea, masterSegmentOrder, slaveSegmentOrder, integrationPoints)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public CellType CellType { get; } = CellType.Unknown;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		private double[] NodalXUpdated()
		{
			var xUpd = new double[2 * Nodes.Count];
			for (var i = 0; i < Nodes.Count; i++)
			{
				xUpd[i * 2] = Nodes[i].X + DisplacementVector[i * 2];
				xUpd[i * 2 + 1] = Nodes[i].Y + DisplacementVector[i * 2 + 1];
			}
			return xUpd;
		}
		private double[] NodalXInitial()
		{
			var x0 = new double[2 * Nodes.Count];
			for (var i = 0; i < Nodes.Count; i++)
			{
				x0[i * 2] = Nodes[i].X;
				x0[i * 2 + 1] = Nodes[i].Y;
			}
			return x0;
		}
		private Tuple<double[,], double[,], double[,], double[,], double[,]> CalculatePositionMatrix(double ksi1, double ksi2)
		{
			if (SlaveSegmentOrder == 1 && MasterSegmentOrder == 1)
			{
				var N1 = 1.0 / 2.0 * (1.0 - ksi1);
				var N2 = 1.0 / 2.0 * (1.0 + ksi1);
				var N3 = 1.0 / 2.0 * (1.0 - ksi2);
				var N4 = 1.0 / 2.0 * (1.0 + ksi2);
				var dN11 = -1.0 / 2.0;
				var dN21 = 1.0 / 2.0;
				var dN32 = -1.0 / 2.0;
				var dN42 = 1.0 / 2.0;
				var aMatrix = new double[,]
				{
					{ -N1, 0.0, -N2, 0.0, N3, 0.0, N4, 0.0 },
					{ 0.0, -N1, 0.0, -N2, 0.0, N3, 0.0, N4 }
				};
				var da1Matrix = new double[,]
				{
					{ -dN11, 0.0, -dN21, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11, 0.0, -dN21, 0.0, 0.0, 0.0, 0.0 }
				};
				var da11Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
				};
				var da2Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, dN32, 0.0, dN42, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, dN32, 0.0, dN42 }
				};
				var da22Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
				};
				return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
			}
			else if (MasterSegmentOrder == 1 && SlaveSegmentOrder == 2)
			{
				var N1 = 1.0 / 2.0 * (1.0 - ksi1);
				var N2 = 1.0 / 2.0 * (1.0 + ksi1);
				var N3 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) - ksi2);
				var N4 = 1.0 - Math.Pow(ksi2, 2.0);
				var N5 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) + ksi2);
				var dN11 = -1.0 / 2.0;
				var dN21 = 1.0 / 2.0;
				var dN32 = ksi2 - 0.5;
				var dN42 = -2.0 * ksi2;
				var dN52 = ksi2 + 0.5;
				var dN322 = 1.0;
				var dN422 = -2.0;
				var dN522 = 1.0;
				var aMatrix = new double[,]
				{
					{ -N1, 0.0, -N2, 0.0, N3, 0.0 , N4, 0.0, N5, 0.0 },
					{ 0.0, -N1, 0.0, -N2, 0.0, N3, 0.0, N4, 0.0, N5 }
				};
				var da1Matrix = new double[,]
				{
					{ -dN11, 0.0, -dN21, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11, 0.0, -dN21, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 }
				};
				var da11Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 }
				};
				var da2Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, dN32, 0.0 , dN42, 0.0, dN52, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, dN32 , 0.0, dN42, 0.0, dN52 }
				};
				var da22Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, dN322, 0.0 , dN422, 0.0, dN522, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, dN322 , 0.0, dN422, 0.0, dN522 }
				};
				return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
			}
			else if (MasterSegmentOrder == 2 && SlaveSegmentOrder == 1)
			{

				var N1 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) - ksi1);
				var N2 = 1.0 - Math.Pow(ksi1, 2.0);
				var N3 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) + ksi1);
				var N4 = 1.0 / 2.0 * (1.0 - ksi2);
				var N5 = 1.0 / 2.0 * (1.0 + ksi2);
				var dN11 = ksi1 - 0.5;
				var dN21 = -2.0 * ksi1;
				var dN31 = ksi1 + 0.5;
				var dN42 = -1.0 / 2.0;
				var dN52 = 1.0 / 2.0;
				var dN111 = 1.0;
				var dN211 = -2.0;
				var dN311 = 1.0;
				var aMatrix = new double[,]
				{
					{ -N1 ,0.0 ,-N2 ,0.0 ,-N3 , 0.0 , N4, 0.0, N5, 0.0 },
					{0.0, -N1 , 0.0 ,-N2, 0.0, -N3, 0.0, N4, 0.0, N5 }
				};
				var da1Matrix = new double[,]
				{
					{ -dN11 ,0.0 ,-dN21 ,0.0 ,-dN31 , 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{0.0, -dN11 , 0.0 ,-dN21, 0.0, -dN31, 0.0, 0.0, 0.0, 0.0 }
				};
				var da11Matrix = new double[,]
				{
					{ -dN111 ,0.0 ,-dN211 ,0.0 ,-dN311 , 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{0.0, -dN111 , 0.0 ,-dN211, 0.0, -dN311, 0.0, 0.0, 0.0, 0.0 }
				};
				var da2Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0 ,0.0 ,0.0 , 0.0 , dN42, 0.0, dN52, 0.0 },
					{0.0, 0.0 , 0.0 ,0.0, 0.0, 0.0, 0.0, dN42, 0.0, dN52 }
				};
				var da22Matrix = new double[,]
				{
					{ 0.0, 0.0 ,0.0 ,0.0 ,0.0 , 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{0.0, 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
				};
				return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
			}
			else
			{
				var N1 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) - ksi1);
				var N2 = 1.0 - Math.Pow(ksi1, 2.0);
				var N3 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) + ksi1);
				var N4 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) - ksi2);
				var N5 = 1.0 - Math.Pow(ksi2, 2.0);
				var N6 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) + ksi2);
				var dN11 = ksi1 - 0.5;
				var dN21 = -2.0 * ksi1;
				var dN31 = ksi1 + 0.5;
				var dN42 = ksi2 - 0.5;
				var dN52 = -2.0 * ksi2;
				var dN62 = ksi2 + 0.5;
				var dN111 = 1.0;
				var dN211 = -2.0;
				var dN311 = 1.0;
				var dN422 = 1.0;
				var dN522 = -2.0;
				var dN622 = 1.0;
				var aMatrix = new double[,]
				{
					{ -N1 ,0.0 ,-N2 ,0.0 ,-N3 , 0.0 , N4, 0.0, N5, 0.0, N6, 0.0 },
					{0.0, -N1 , 0.0 ,-N2, 0.0, -N3, 0.0, N4, 0.0, N5, 0.0, N6 }
				};
				var da1Matrix = new double[,]
				{
					{ -dN11 ,0.0 ,-dN21 ,0.0 ,-dN31, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{0.0, -dN11 , 0.0 ,-dN21, 0.0, -dN31, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
				};
				var da11Matrix = new double[,]
				{
					{ -dN111 ,0.0 ,-dN211 ,0.0 ,-dN311, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{0.0, -dN111 , 0.0 ,-dN211, 0.0, -dN311, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
				};
				var da2Matrix = new double[,]
				{
					{ 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, dN42, 0.0, dN52, 0.0, dN62, 0.0 },
					{ 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, dN42, 0.0, dN52, 0.0, dN62 }
				};
				var da22Matrix = new double[,]
				{
					{ 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, dN422, 0.0, dN522, 0.0, dN622, 0.0 },
					{ 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, dN422, 0.0, dN522, 0.0, dN622 }
				};
				return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
			}
		}

		private Tuple<double[], double, double[], double[], double, double> MasterSegmentGeometry(double[,] daMatrix, double[,] da2Matrix)
		{
			var xupd = NodalXUpdated().Scale(-1.0);
			var surfaceVector = Matrix.CreateFromArray(daMatrix).Multiply(xupd);
			var surfaceVectorDerivative = Matrix.CreateFromArray(da2Matrix).Multiply(xupd);
			var detm = surfaceVector.DotProduct(surfaceVector);
			var m11 = 1.0 / detm;
			var vector = new double[2];
			var scalarCoef = new double();
			var scalarCoef2 = Math.Pow(m11, 2.0);
			var tangentVector = new double[2];
			if (MasterSegmentOrder == 1)
			{
				var Xm1 = Nodes[0].X + DisplacementVector[0];
				var Ym1 = Nodes[0].Y + DisplacementVector[1];
				var Xm2 = Nodes[1].X + DisplacementVector[2];
				var Ym2 = Nodes[1].Y + DisplacementVector[3];
				vector[0] = Ym2 - Ym1;
				vector[1] = Xm1 - Xm2;
				scalarCoef = -1.0 / (2.0 * Math.Sqrt(detm));
			}
			else
			{
				vector[0] = -surfaceVector[1];
				vector[1] = surfaceVector[0];
				scalarCoef = 1.0 / (Math.Sqrt(detm));
				tangentVector = surfaceVector.Scale(scalarCoef);
			}
			var normalUnitVec = vector.Scale(scalarCoef);
			var curvatureTensor = scalarCoef2 * surfaceVectorDerivative.DotProduct(normalUnitVec);
			return new Tuple<double[], double, double[], double[], double, double>(surfaceVector, m11, normalUnitVec, tangentVector, curvatureTensor, detm);
		}

		private Tuple<double[], double, double[], double[], double> SlaveSegmentGeometry(double[,] daMatrix, double[,] da2Matrix)
		{
			var xupd = NodalXUpdated();
			var surfaceVector = Matrix.CreateFromArray(daMatrix).Multiply(xupd);
			var surfaceVectorDerivative = Matrix.CreateFromArray(da2Matrix).Multiply(xupd);
			var detm = surfaceVector.DotProduct(surfaceVector);
			var m11 = 1.0 / detm;
			var vector = new double[] { -surfaceVector[1], surfaceVector[0] };
			var scalarCoef = 1.0 / (Math.Sqrt(detm));
			var normalUnitVec = vector.Scale(scalarCoef);
			var tangentVector = surfaceVector.Scale(scalarCoef);
			var scalarCoef2 = Math.Pow(m11, 2.0);
			var curvatureTensor = scalarCoef2 * surfaceVectorDerivative.DotProduct(normalUnitVec);
			return new Tuple<double[], double, double[], double[], double>(surfaceVector, detm, normalUnitVec, tangentVector, curvatureTensor);
		}

		private double CalculatePenetration(double[,] aMatrix, double[] n)
		{
			var AT = Matrix.CreateFromArray(aMatrix).Transpose();
			var AT_n = AT.Multiply(n);
			var xupd = NodalXUpdated();
			var normalGap = xupd.DotProduct(AT_n);
			return normalGap;
		}

		private double CalculateDeltaKsi(double[] masterSlaveRelativeVector, double[] surfaceVector, double[] surfaceVectorDerivative)
		{
			var scalar1 = surfaceVector.DotProduct(masterSlaveRelativeVector);
			var scalar2 = surfaceVectorDerivative.DotProduct(masterSlaveRelativeVector) - surfaceVector.DotProduct(surfaceVector);
			var deltaKsi = -scalar1 / scalar2;
			return deltaKsi;
		}

		private double Project(double ksi1Initial, double ksi2)
		{
			if (MasterSegmentOrder == 1)
			{
				var aMatrix = CalculatePositionMatrix(ksi1Initial, ksi2).Item1;
				var m = SlaveSegmentOrder + 1;
				var slaveNMatrix = new double[2, 2 * m];
				var xUpdated = NodalXUpdated();
				var list = new List<double>();
				for (var i = 4; i < list.Count; i++)
				{
					list.Add(xUpdated[i]);
				}
				var x = list.ToArray();
				for (var i = 0; i <= 1; i++)
				{
					if (i == 0)
					{
						var countCols = 0;
						for (var j = 4; j < aMatrix.GetLength(1) - 1; j += 2)
						{
							slaveNMatrix[i, countCols] = aMatrix[i, j];
							countCols += 2;
						}
					}
					else
					{
						var countCols = 1;
						for (var j = 5; j < aMatrix.GetLength(1); j += 2)
						{
							slaveNMatrix[i, countCols] = aMatrix[i, j];
							countCols += 2;
						}
					}
				}
				var slavePositionVector = Matrix.CreateFromArray(slaveNMatrix).Multiply(x);
				var xM1 = xUpdated[0];
				var yM1 = xUpdated[1];
				var xM2 = xUpdated[2];
				var yM2 = xUpdated[3];
				var xS = slavePositionVector[0];
				var yS = slavePositionVector[1];
				var ksi = (2 * (xS * (xM2 - xM1) + yS * (yM2 - yM1)) - Math.Pow(xM2, 2) - Math.Pow(yM2, 2) + Math.Pow(xM1, 2) + Math.Pow(yM1, 2)) / (Math.Pow(xM2 - xM1, 2) + Math.Pow(yM2 - yM1, 2));
				return ksi;
			}
			else
			{
				var maxIterations = 1000;
				var tol = Math.Pow(10.0, -6.0);
				var deltaKsi = 0.0;
				var ksi = ksi1Initial;
				var xUpdated = NodalXUpdated();
				for (var i = 1; i <= maxIterations; i++)
				{
					var aMatrices = CalculatePositionMatrix(ksi, ksi2);
					var masterSlaveRelativeVector = Matrix.CreateFromArray(aMatrices.Item1).Multiply(xUpdated);
					var surfaceVector = Matrix.CreateFromArray(aMatrices.Item2).Multiply(xUpdated).Scale(-1.0);
					var surfaceVectorDerivative = Matrix.CreateFromArray(aMatrices.Item3).Multiply(xUpdated).Scale(-1.0);
					deltaKsi = CalculateDeltaKsi(masterSlaveRelativeVector, surfaceVector, surfaceVectorDerivative);
					ksi += deltaKsi;
					if (Math.Abs(deltaKsi) <= tol)
					{
						break;
					}
				}
				if (Math.Abs(deltaKsi) > tol)
				{
					throw new Exception("CPP not found in current iterations");
				}
				else
				{
					return ksi;
				}
			}
		}
		private double ProjectInitialConfiguration(double ksi1Initial, double ksi2)
		{
			var maxIterations = 1000;
			var tol = Math.Pow(10.0, -6.0);
			var deltaKsi = 0.0;
			var ksi = ksi1Initial;
			var xInitial = NodalXInitial();
			for (var i = 1; i <= maxIterations; i++)
			{
				var aMatrices = CalculatePositionMatrix(ksi, ksi2);
				var masterSlaveRelativeVector = Matrix.CreateFromArray(aMatrices.Item1).Multiply(xInitial);
				var surfaceVector = (Matrix.CreateFromArray(aMatrices.Item2).Multiply(xInitial)).Scale(-1.0);
				var surfaceVectorDerivative = (Matrix.CreateFromArray(aMatrices.Item3).Multiply(xInitial)).Scale(-1.0);
				deltaKsi = CalculateDeltaKsi(masterSlaveRelativeVector, surfaceVector, surfaceVectorDerivative);
				ksi += deltaKsi;
				if (Math.Abs(deltaKsi) <= tol)
				{
					break;
				}
			}
			if (Math.Abs(deltaKsi) > tol)
			{
				throw new Exception("CPP not found in current iterations");
			}
			else
			{
				return ksi;
			}
		}
		public void InitializeTangentialProperties()
		{
			var gPointsArray = GaussPoints().Item1;
			for (var i = 0; i < IntegrationPoints; i++)
			{
				var ksi2 = gPointsArray[i];
				var ksi11 = ProjectInitialConfiguration(0.0, ksi2);
				IntegrationPointsStickingPoints[i + 1] = ksi11;
			}
		}
		public void UpdateTangentialProperties()
		{
			var gPointsArray = GaussPoints().Item1;
			for (var i = 0; i < IntegrationPoints; i++)
			{
				var ksi2 = gPointsArray[i];
				var ksi11 = Project(0.0, ksi2);
				var deltaKsi1 = ksi11 - IntegrationPointsStickingPoints[i + 1];
				var positionMatrices = CalculatePositionMatrix(ksi11, ksi2);
				var da1Matrix = positionMatrices.Item2;
				var da11Matrix = positionMatrices.Item3;
				var detm = MasterSegmentGeometry(da1Matrix, da11Matrix).Item6;
				var tangentialTraction = IntegrationPointsTangentialTractions[i + 1] - detm * PenaltyFactorTangential * deltaKsi1;
				IntegrationPointsStickingPoints[i + 1] = ksi11;
				IntegrationPointsTangentialTractions[i + 1] = tangentialTraction;
			}
		}

		private Tuple<double[], double[]> GaussPoints()
		{
			var iP = IntegrationPoints;
			var gaussPoints = new double[iP];
			var gaussWeights = new double[iP];
			if (iP == 2)
			{
				gaussPoints[0] = -1.0 / Math.Sqrt(3);
				gaussPoints[1] = 1.0 / Math.Sqrt(3);
				gaussWeights[0] = 1.0;
				gaussWeights[1] = 1.0;
			}
			else if (iP == 3)
			{
				gaussPoints[0] = -0.77459;
				gaussPoints[1] = 0.0;
				gaussPoints[2] = 0.77459;
				gaussWeights[0] = 0.55555;
				gaussWeights[1] = 0.88888;
				gaussWeights[2] = 0.55555;
			}
			else if (iP == 4)
			{
				gaussPoints[0] = -0.86113;
				gaussPoints[1] = -0.33998;
				gaussPoints[2] = 0.33998;
				gaussPoints[3] = 0.86113;
				gaussWeights[0] = 0.34785;
				gaussWeights[1] = 0.65214;
				gaussWeights[2] = 0.65214;
				gaussWeights[3] = 0.34785;
			}
			else if (iP == 5)
			{
				gaussPoints[0] = -0.90617;
				gaussPoints[1] = -0.53846;
				gaussPoints[2] = 0.0;
				gaussPoints[3] = 0.53846;
				gaussPoints[4] = 0.90617;
				gaussWeights[0] = 0.23692;
				gaussWeights[1] = 0.47862;
				gaussWeights[2] = 0.56888;
				gaussWeights[3] = 0.47862;
				gaussWeights[4] = 0.23692;
			}
			else if (iP == 6)
			{
				gaussPoints[0] = -0.93246;
				gaussPoints[1] = -0.66120;
				gaussPoints[2] = -0.23861;
				gaussPoints[3] = 0.23861;
				gaussPoints[4] = 0.66120;
				gaussPoints[5] = 0.93246;
				gaussWeights[0] = 0.17132;
				gaussWeights[1] = 0.36076;
				gaussWeights[2] = 0.46791;
				gaussWeights[3] = 0.46791;
				gaussWeights[4] = 0.36076;
				gaussWeights[5] = 0.17132;
			}
			else if (iP == 7)
			{
				gaussPoints[0] = -0.94910;
				gaussPoints[1] = -0.74153;
				gaussPoints[2] = -0.40584;
				gaussPoints[3] = 0.0;
				gaussPoints[4] = 0.40584;
				gaussPoints[5] = 0.74153;
				gaussPoints[6] = 0.94910;
				gaussWeights[0] = 0.12948;
				gaussWeights[1] = 0.27970;
				gaussWeights[2] = 0.38183;
				gaussWeights[3] = 0.41795;
				gaussWeights[4] = 0.38183;
				gaussWeights[5] = 0.27970;
				gaussWeights[6] = 0.12948;
			}
			else if (iP == 8)
			{
				gaussPoints[0] = -0.96028;
				gaussPoints[1] = -0.79666;
				gaussPoints[2] = -0.52553;
				gaussPoints[3] = -0.18343;
				gaussPoints[4] = 0.18343;
				gaussPoints[5] = 0.52553;
				gaussPoints[6] = 0.79666;
				gaussPoints[7] = 0.96028;
				gaussWeights[0] = 0.10122;
				gaussWeights[1] = 0.22238;
				gaussWeights[2] = 0.31370;
				gaussWeights[3] = 0.36268;
				gaussWeights[4] = 0.36268;
				gaussWeights[5] = 0.31370;
				gaussWeights[6] = 0.22238;
				gaussWeights[7] = 0.10122;
			}
			else if (iP == 9)
			{
				gaussPoints[0] = -0.96816;
				gaussPoints[1] = -0.83603;
				gaussPoints[2] = -0.61337;
				gaussPoints[3] = -0.32425;
				gaussPoints[4] = 0.0;
				gaussPoints[5] = 0.32425;
				gaussPoints[6] = 0.61337;
				gaussPoints[7] = 0.83603;
				gaussPoints[8] = 0.96816;
				gaussWeights[0] = 0.08127;
				gaussWeights[1] = 0.18064;
				gaussWeights[2] = 0.26061;
				gaussWeights[3] = 0.31234;
				gaussWeights[4] = 0.33023;
				gaussWeights[5] = 0.31234;
				gaussWeights[6] = 0.26061;
				gaussWeights[7] = 0.18064;
				gaussWeights[8] = 0.08127;
			}
			else if (iP == 10)
			{
				gaussPoints[0] = -0.97390;
				gaussPoints[1] = -0.86506;
				gaussPoints[2] = -0.67940;
				gaussPoints[3] = -0.43339;
				gaussPoints[4] = -0.14887;
				gaussPoints[5] = 0.14887;
				gaussPoints[6] = 0.43339;
				gaussPoints[7] = 0.67940;
				gaussPoints[8] = 0.86506;
				gaussPoints[9] = 0.97390;
				gaussWeights[0] = 0.06667;
				gaussWeights[1] = 0.14945;
				gaussWeights[2] = 0.21908;
				gaussWeights[3] = 0.26926;
				gaussWeights[4] = 0.29552;
				gaussWeights[5] = 0.29552;
				gaussWeights[6] = 0.26926;
				gaussWeights[7] = 0.21908;
				gaussWeights[8] = 0.14945;
				gaussWeights[9] = 0.06667;
			}
			return new Tuple<double[], double[]>(gaussPoints, gaussWeights);
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		private Matrix CalculateMainStiffnessPart(Matrix A, double[] n)
		{
			var nxn = n.TensorProduct(n);
			var nxn_A = nxn * A;
			var AT_nxn_A = A.Transpose() * nxn_A;
			var mainStiffnessMatrix = AT_nxn_A.Scale(PenaltyFactorNormal);
			return mainStiffnessMatrix;
		}

		private Matrix CalculateRotationalStiffnessPart(double[,] A, double[,] dA, double[] n, double ksi3, double m11, double[] dRho)
		{
			var coef = PenaltyFactorNormal * ksi3 * m11;
			var n_x_dRho = n.TensorProduct(dRho);
			var dRho_x_n = dRho.TensorProduct(n);
			var firstTerm = Matrix.CreateFromArray(dA).Transpose() * n_x_dRho * Matrix.CreateFromArray(A);
			var secondTerm = Matrix.CreateFromArray(A).Transpose() * dRho_x_n * Matrix.CreateFromArray(dA);
			var rotationalPart = (firstTerm + secondTerm).Scale(coef);
			return rotationalPart;
		}

		private Matrix CalculateCurvatureStiffnessPart(double[,] A, double ksi3, double m11, double[] dRho, double h11)
		{
			var coef = PenaltyFactorNormal * ksi3 * m11 * h11;
			var dRho_x_dRho = dRho.TensorProduct(dRho);
			var mat = Matrix.CreateFromArray(A).Transpose() * dRho_x_dRho * Matrix.CreateFromArray(A);
			var curvaturePart = mat.Scale(coef);
			return curvaturePart;
		}

		private Matrix CalculateTangentialStiffnessPartForSticking(Matrix A, Matrix dA, double Tangtr, double[] normalVector, double[] tangentVector, double detm, double m11, double h11)
		{
			var coef2 = Tangtr * m11;
			var coef3 = Tangtr * h11 * Math.Pow(detm, 0.5);
			var t_x_t = tangentVector.TensorProduct(tangentVector);
			var t_x_n = tangentVector.TensorProduct(normalVector);
			var n_x_t = normalVector.TensorProduct(tangentVector);
			var tangentialPart1 = (A.Transpose() * t_x_t * A).Scale(PenaltyFactorTangential);
			var tangentialPart2 = (dA.Transpose() * t_x_t * A).Scale(coef2);
			var tangentialPart3 = (A.Transpose() * (t_x_n + n_x_t) * A).Scale(coef3);
			var fullTangentialPart = tangentialPart1.Scale(-1d) + tangentialPart2 + tangentialPart2.Transpose() + tangentialPart3.Transpose();
			return fullTangentialPart;
		}
		private Matrix CalculateTangentialStiffnessPartForSliding(Matrix A, Matrix dA, double Tangtr, double[] normalVector, double[] tangentVector, double detm, double m11, double h11, double ksi3)
		{
			var coef1 = SlidingCoefficient * PenaltyFactorNormal * (Tangtr / Math.Abs(Tangtr));
			var coef2 = SlidingCoefficient * PenaltyFactorNormal * (Tangtr / Math.Abs(Tangtr)) * Math.Abs(ksi3) * Math.Pow(m11, 0.5);
			var coef3 = SlidingCoefficient * PenaltyFactorNormal * (Tangtr / Math.Abs(Tangtr)) * Math.Abs(ksi3) * h11 * detm;
			var t_x_t = tangentVector.TensorProduct(tangentVector);
			var t_x_n = tangentVector.TensorProduct(normalVector);
			var n_x_t = normalVector.TensorProduct(tangentVector);
			var tangentialPart1 = (A.Transpose() * t_x_n * A).Scale(coef1);
			var tangentialPart2 = (dA.Transpose() * t_x_t * A).Scale(coef2);
			var tangentialPart3 = (A.Transpose() * (t_x_n + n_x_t) * A).Scale(coef3);
			var fullTangentialPart = tangentialPart1.Scale(-1d) + tangentialPart2 + tangentialPart2.Transpose() + tangentialPart3.Transpose();
			return fullTangentialPart;
		}

		public IMatrix StiffnessMatrix()
		{
			var globalStifnessMatrix = Matrix.CreateFromArray(new double[2 * Nodes.Count, 2 * Nodes.Count]);
			var gPArray = GaussPoints().Item1;
			var gWArray = GaussPoints().Item2;
			for (var i = 0; i < IntegrationPoints; i++)
			{
				var ksi2 = gPArray[i];
				var gW = gWArray[i];
				var ksi1 = Project(0.0, ksi2);
				if (Math.Abs(ksi1) <= 1.05)
				{
					var positionMatrices = CalculatePositionMatrix(ksi1, ksi2);
					var aMatrix = positionMatrices.Item1;
					var da1Matrix = positionMatrices.Item2;
					var da11Matrix = positionMatrices.Item3;
					var masterSurfaceCharacteristics = MasterSegmentGeometry(da1Matrix, da11Matrix);
					var dRho = masterSurfaceCharacteristics.Item1;
					var m11 = masterSurfaceCharacteristics.Item2;
					var n = masterSurfaceCharacteristics.Item3;
					var t = masterSurfaceCharacteristics.Item4;
					var h11 = masterSurfaceCharacteristics.Item5;
					var detm = masterSurfaceCharacteristics.Item6;
					var ksi3 = CalculatePenetration(aMatrix, n);
					if (ksi3 <= 0)
					{
						var slaveMetricTensor = SlaveSegmentGeometry(positionMatrices.Item4, positionMatrices.Item5).Item2;
						var mainPart = CalculateMainStiffnessPart(Matrix.CreateFromArray(aMatrix), n);
						var rotationalPart = CalculateRotationalStiffnessPart(aMatrix, da1Matrix, n, ksi3, m11, dRho);
						var curvaturePart = CalculateCurvatureStiffnessPart(aMatrix, ksi3, m11, dRho, h11);
						var scalar = Math.Pow(slaveMetricTensor, 0.5) * gW;
						var stifnessMatrixNormalPartIntegrationPointContribution = (mainPart + rotationalPart + curvaturePart).Scale(scalar);
						globalStifnessMatrix += stifnessMatrixNormalPartIntegrationPointContribution;
						var deltaKsi1 = ksi1 - IntegrationPointsStickingPoints[i + 1];
						var trialTangentialTraction = IntegrationPointsTangentialTractions[i + 1] - detm * PenaltyFactorTangential * deltaKsi1;
						var phiTr = Math.Pow(trialTangentialTraction * trialTangentialTraction * m11, 0.5) - StickingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3);
						if (phiTr <= 0)
						{
							//Sticking
							var StifnessMatrixTangentialPart = CalculateTangentialStiffnessPartForSticking(Matrix.CreateFromArray(aMatrix), Matrix.CreateFromArray(da1Matrix),
								trialTangentialTraction, n, t, detm, m11, h11).Scale(scalar);
							globalStifnessMatrix += StifnessMatrixTangentialPart;
						}
						else
						{
							//Sliding
							var StifnessMatrixTangentialPart = CalculateTangentialStiffnessPartForSliding(Matrix.CreateFromArray(aMatrix), Matrix.CreateFromArray(da1Matrix),
								trialTangentialTraction, n, t, detm, m11, h11, ksi3).Scale(scalar);
							globalStifnessMatrix += StifnessMatrixTangentialPart;
						}
					}
				}
			}
			return dofEnumerator.GetTransformedMatrix(globalStifnessMatrix);
		}

		public IMatrix PhysicsMatrix() => StiffnessMatrix();

		public IMatrix MassMatrix()
		{
			double[,] massMatrix = new double[2 * Nodes.Count, 2 * Nodes.Count];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(massMatrix));
		}

		public IMatrix DampingMatrix()
		{
			double[,] dampingMatrix = new double[2 * Nodes.Count, 2 * Nodes.Count];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(dampingMatrix));
		}

		private double[] CreateInternalGlobalForcesVector()
		{
			var internalGlobalForcesVector = new double[2 * Nodes.Count];
			var gPArray = GaussPoints().Item1;
			var gWArray = GaussPoints().Item2;
			for (var i = 0; i < IntegrationPoints; i++)
			{
				var ksi2 = gPArray[i];
				var gW = gWArray[i];
				var ksi1 = Project(0.0, ksi2);
				if (Math.Abs(ksi1) <= 1.05)
				{
					var positionMatrices = CalculatePositionMatrix(ksi1, ksi2);
					var aMatrix = positionMatrices.Item1;
					var da1Matrix = positionMatrices.Item2;
					var da11Matrix = positionMatrices.Item3;
					var surfaceCharacteristics = MasterSegmentGeometry(da1Matrix, da11Matrix);
					var m11 = surfaceCharacteristics.Item2;
					var n = surfaceCharacteristics.Item3;
					var t = surfaceCharacteristics.Item4;
					var detm = surfaceCharacteristics.Item6;
					var ksi3 = CalculatePenetration(aMatrix, n);
					if (ksi3 <= 0)
					{
						var slaveMetricTensor = SlaveSegmentGeometry(positionMatrices.Item4, positionMatrices.Item5).Item2;
						var scalar = Math.Pow(slaveMetricTensor, 0.5) * gW;
						var AT = Matrix.CreateFromArray(aMatrix).Transpose();
						var AT_n = AT.Multiply(n);
						var internalForcesVectorGaussPointNormalPart = AT_n.Scale(PenaltyFactorNormal * ksi3 * scalar);
						internalGlobalForcesVector = internalGlobalForcesVector.Add(internalForcesVectorGaussPointNormalPart);
						var deltaKsi1 = ksi1 - IntegrationPointsStickingPoints[i + 1];
						var trialTangentialTraction = IntegrationPointsTangentialTractions[i + 1] - detm * PenaltyFactorTangential * deltaKsi1;
						var phiTr = Math.Pow(trialTangentialTraction * trialTangentialTraction * m11, 0.5) - StickingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3);
						if (phiTr <= 0)
						{
							//Sticking
							var AT_t = AT.Multiply(t);
							var internalTangentialForcesVectorGaussPoint = AT_t.Scale(trialTangentialTraction * Math.Pow(m11, 0.5) * scalar);
							internalGlobalForcesVector = internalGlobalForcesVector.Add(internalTangentialForcesVectorGaussPoint);
						}
						else
						{
							//Sliding
							var T1 = (trialTangentialTraction / Math.Abs(trialTangentialTraction)) * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3) * Math.Pow(detm, 0.5);
							var AT_t = AT.Multiply(t);
							var internalTangentialForcesVectorGaussPoint = AT_t.Scale(T1 * Math.Pow(m11, 0.5) * scalar);
							internalGlobalForcesVector = internalGlobalForcesVector.Add(internalTangentialForcesVectorGaussPoint);
						}
					}
				}
			}
			return internalGlobalForcesVector;
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
