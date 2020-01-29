using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Integration.Curves.ErrorReporting;

namespace Integration.Curves
{
    public class MonotonicallyIncreasing : PairedData
    {
        #region Notes
        #endregion
        #region Fields
        #endregion
        #region Properties
        #endregion
        #region Constructors
        public MonotonicallyIncreasing(List<double> XValues, List<double> YValues)
        {
            X = XValues.ToArray();
            Y = YValues.ToArray();
        }


        #endregion
        #region Voids
        #endregion
        #region Functions
        public override double GetXFromY(double y)
        {
            throw new NotImplementedException();
        }

        public override double GetYFromX(double x)
        {
            throw new NotImplementedException();
        }

        public override ErrorReport Validate()
        {
            ErrorReporting.ErrorReport report = new ErrorReport();
            if (X.Count() != Y.Count()) { report.AddError(new CurveError("Your X values and Y values are not the same length.", 0, 0)); }
            int increasingCount = 0;
            for(int i = 1; i< Math.Min(X.Count(), Y.Count()); i++)
            {
                if (Y[i] < Y[i - 1]) { report.AddError(new CurveError("Input Y values are not monotonically increasing",i,1)); }
                if (X[i] < X[i - 1]) { report.AddError(new CurveError("Input Y values are not monotonically increasing", i, 1)); increasingCount += 1; }
                if(increasingCount>0 && increasingCount != Math.Min(Y.Count(), X.Count()) - 2) { report.AddError(new CurveError("The input arrays do not pass the vertical line test", 0, 0)); }
            }
            return report;
        }
        #endregion
    }
}
