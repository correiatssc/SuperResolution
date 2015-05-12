using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SuperResolution
{
    public class ExcecaoEsperada : Exception
    {
        public ExcecaoEsperada(String mensagem)
            : base(mensagem)
        { }
    }
}
