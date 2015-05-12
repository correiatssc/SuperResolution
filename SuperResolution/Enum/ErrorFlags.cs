using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SuperResolution.Enum
{
    public enum ErrorFlags
    {
        CANT_WRITE, 
        DECI_FCTR, 
        INSUF_MEM, 
        FILE_ERR, 
        FILE_SEEK, 
        MAX_FREQ, 
        MIN_FREQ,
        MISUSE, 
        NOISE_FLOOR, 
        SAMPLE_FREQ, 
        SFT_OOR, 
        THR_DH, 
        THR_HIGH, 
        THR_MAX_RTO,
        THR_MIN
    }
}
