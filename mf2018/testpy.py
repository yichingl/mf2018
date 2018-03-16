# test matlab engine first by doing a local call, then by passing engine to helper function

import matlab.engine
import testpy2
eng = matlab.engine.start_matlab();

print("testing hello");
ret = testpy2.hello_fn();
print(ret);

print("testing local call");
ret = eng.test(11);
print(ret);

print("testing by passing engine");
ret = testpy2.test_fn(eng,12);
print(ret);