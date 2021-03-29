import json,click,math
import matplotlib.pyplot as plt

def plot(series,title):
    fig = plt.figure()
    fig.suptitle(title)
    ax = fig.subplots()
    for s in series:
        x = [ xy[0] for xy in series[s]]
        y = [ xy[1] for xy in series[s]]
        #ax.semilogx(x,y,'-o',label=s)
        ax.loglog(x,y,'-o',label=s)
    ax.set_xlabel('size')
    ax.set_ylabel('speed (mflops)')
    ax.legend()    
    plt.show()    

@click.command()
@click.argument('fname')
@click.argument('title')
def main(fname,title):
    db = json.load(open(fname,'r'))
    series = { 
        'fftx::InPlace' : [], 
        'fftx::DivideAndConquer' : [], 
        'fftx::Iterative' : [] ,
        'FFTW' : [] }
    def is_fftx(name):
        return name=='InPlace' or name == 'DivideAndConquer' or name == 'Iterative'
    def is_fftw(name):
        return name=='FFTW'
        
    for b in db['benchmarks']:
        try:
            name,val = b['name'][6:].split('/')
            val = int(val)
            time = float(b['real_time'])/1e9 # seconds
            mflops = 5*val*math.log(val)/math.log(2)/time/1e6
            if is_fftx(name):
                name = 'fftx::' + name
                series[name].append((val,mflops))
            if is_fftw(name):
                series[name].append((val,mflops))
        except:
            pass
    plot(series,title)    
    

if __name__=="__main__":
    main()
