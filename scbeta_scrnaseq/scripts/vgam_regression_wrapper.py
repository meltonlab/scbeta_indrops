import os,sys
import loompy

if __name__=="__main__":
    sys.path.append(os.path.realpath(os.path.join(sys.argv[0], '..', '..', '..')))
    print(os.path.realpath(os.path.join(sys.argv[0], '..', '..', '..')))
    import scbeta_scrnaseq.pseudotime

    ds =  loompy.connect(sys.argv[1])
    scbeta_scrnaseq.pseudotime.vgam_regression(ds)
    ds.close()
