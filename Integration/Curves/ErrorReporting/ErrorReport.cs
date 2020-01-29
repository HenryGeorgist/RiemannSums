using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Integration.Curves.ErrorReporting
{
    public class ErrorReport
    {

        #region Properties
        public List<CurveError> Errors { get; set; }
        #endregion
        #region Constructors
        public ErrorReport(List<CurveError> errs)
        {
            Errors = errs;
        }
        public ErrorReport()
        {
            Errors = new List<CurveError>();
        }
        #endregion
        #region Voids
        public void AddError(CurveError err)
        {
            Errors.Add(err);
        }
        #endregion
        #region Functions
        #endregion
    }
}
