### Hints and Builtins research for emulated modular arithmetic in Cairo. 

setup for research : 
```
python3 -m venv venv
echo 'export PYTHONPATH="$PWD:$PYTHONPATH"' >> venv/bin/activate
source venv/bin/activate
```

next times : 
```
source venv/bin/activate
```

#### Current status : 

- [src] Mul Hint written in Rust ok.
- [research] AIR and polynomial constraints for mul (without range checks) ok. 