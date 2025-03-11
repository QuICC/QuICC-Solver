import io
from contextlib import redirect_stdout
import re

def cppify(func):
    sout = io.StringIO()
    with redirect_stdout(sout):
       func()
    cpp_out = sout.getvalue()
    print(cpp_out)

    cpp_out = re.sub(': \t', ': return ', cpp_out)
    cpp_out = re.sub('\t', '', cpp_out)
    cpp_out = re.sub("([-]?\\d+):", '[\\g<1>]:', cpp_out)
    cpp_out = re.sub("\\*\\*(\\d+)", '<\\g<1>>()', cpp_out)
    cpp_out = re.sub("([ \\(-]\\d+(?!]))", '\\1.0_mp', cpp_out)
    cpp_out = re.sub("([ \\*\\(][abl])([ \\*\\)/])", "\\g<1><1>()\\g<2>", cpp_out)
    cpp_out = re.sub("\\n", ";\\n", cpp_out)
    cpp_out = re.sub("n<2>\\(\\)", "n*n", cpp_out)
    cpp_out = re.sub("n<3>\\(\\)", "n*n*n", cpp_out)
    cpp_out = re.sub("n<4>\\(\\)", "n*n*n*n", cpp_out)
    cpp_out = re.sub("n<5>\\(\\)", "n*n*n*n*n", cpp_out)

    print(cpp_out)
