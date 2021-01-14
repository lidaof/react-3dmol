const resnList = ['GLN','THR','SER','VAL','PRO','LYS','ILE','LEU','ARG','GLY','CYS','ASP','GLU','ASN','TYR','MET','ALA','PHE','TRP','HIS'];


export const g3dParser = function(data) {
    /**
     * the .g3d parser for 3dmol.js input, following .xyz example
     */
    const atoms = []; //array of models, maybe for different haplotypes, or cells/clusters for single cell data
    for (const hap of Object.keys(data)) {
      const hapAtoms = [];
      let bi = 0 // bond index
      for (const chr of Object.keys(data[hap])) {
        for(let i=0; i < data[hap][chr].start.length; i++) {
            bi += 1;
          const atom = {};
          const start = data[hap][chr].start[i];
          atom.x = data[hap][chr].x[i];
          atom.y = data[hap][chr].y[i];
          atom.z = data[hap][chr].z[i];
          atom.resi = bi;
          atom.resn = resnList[i%resnList.length]; 
          atom.atom = 'CA'
          atom.elem = 'Ca'
          atom.chain = 'A';
          atom.chrom = chr;
          atom.serial = bi;
          atom.b = 1.0
        //   atom.atom = `${chr}-${i}`;
        //   atom.resn = `${chr}-${start}`;
          atom.hetflag = false;
          atom.bonds = [];
          atom.bondOrder = [];
          atom.properties = {};
    //       if(i === 0) {s
    //         atom.bonds = [bi, bi+1, bi+2, bi+3];
    //     }else if(i === 1) {
    //       atom.bonds = [bi-2, bi, bi+1, bi+2];
    //   }else if(i === data[hap][chr].start.length-1) {
    //       atom.bonds = [bi-2, bi-3, bi-4, bi-5];
    //     }else if(i === data[hap][chr].start.length-2) {
    //       atom.bonds = [bi, bi-2, bi-3, bi-4];
    //     } else if (i>1 && i<data[hap][chr].start.length-2) {
    //       atom.bonds = [bi-3, bi-2, bi, bi+1];
    //   }
    //   atom.bondOrder = [1,1,1,1]
//       if(i === 0) {
//         atom.bonds = [bi];
//         atom.bondOrder = [1]
//     }else if(i === data[hap][chr].start.length-1) {
//       atom.bonds = [bi-2];
//       atom.bondOrder = [1]
//     }else {
//       atom.bonds = [ bi-2, bi];
//       atom.bondOrder = [1, 1];
//   }
        //   atom.properties.value = Math.floor(Math.random() * 201) - 100; // random value for test color assignment
        //   atom.properties.chrom = chr;
        //   atom.properties.start = start;
        //   atom.properties.radius = 50;
          hapAtoms.push(atom);
        }
      }
      atoms.push(hapAtoms);
    }
    for (let i = 0; i < atoms.length; i++) {
      assignBonds(atoms[i]);
    }
    console.log(atoms)
     return atoms;
  };


const bondTable = {
    H :0.37,                                                                                                                                He:0.32,
    Li:1.34,Be:0.90,                                                                                B :0.82,C :0.77,N :0.75,O :0.73,F :0.71,Ne:0.69,
    Na:1.54,Mg:1.30,                                                                                Al:1.18,Si:1.11,P :1.06,S :1.02,Cl:0.99,Ar:0.97,
    K :1.96,Ca:1.74,Sc:1.44,Ti:1.56,V :1.25,/* Cr */Mn:1.39,Fe:1.25,Co:1.26,Ni:1.21,Cu:1.38,Zn:1.31,Ga:1.26,Ge:1.22,/* As */Se:1.16,Br:1.14,Kr:1.10,
    Rb:2.11,Sr:1.92,Y :1.62,Zr:1.48,Nb:1.37,Mo:1.45,Tc:1.56,Ru:1.26,Rh:1.35,Pd:1.31,Ag:1.53,Cd:1.48,In:1.44,Sn:1.41,Sb:1.38,Te:1.35,I :1.33,Xe:1.30,
    Cs:2.25,Ba:1.98,Lu:1.60,Hf:1.50,Ta:1.38,W :1.46,Re:1.59,Os:1.44,Ir:1.37,Pt:1.28,Au:1.44,Hg:1.49,Tl:1.48,Pb:1.47,Bi:1.46,/* Po *//* At */Rn:1.45,

    // None of the bottom row or any of the Lanthanides have bond lengths
};

const bondLength = function(elem) {
    return bondTable[elem] || 1.6;
};

const areConnected = function(atom1, atom2) {
    var maxsq = bondLength(atom1.elem) + bondLength(atom2.elem);
    maxsq += 0.25;// fudge factor, especially important for md frames, also see 1i3d
    maxsq *= maxsq;

    var xdiff = atom1.x - atom2.x;
    xdiff *= xdiff;
    if (xdiff > maxsq)
        return false;
    var ydiff = atom1.y - atom2.y;
    ydiff *= ydiff;
    if (ydiff > maxsq)
        return false;
    var zdiff = atom1.z - atom2.z;
    zdiff *= zdiff;
    if (zdiff > maxsq)
        return false;

    var distSquared = xdiff + ydiff + zdiff;

    if (isNaN(distSquared))
        return false;
    else if (distSquared < 0.5)
        return false; // maybe duplicate position.
    else if (distSquared > maxsq)
        return false;
    else if(atom1.altLoc !== atom2.altLoc && atom1.altLoc !== ' ' && atom2.altLoc !== ' ')
        return false; // don't connect across alternate locations
    else
        return true;
};
/**
     * @param {AtomSpec[]}
     *            atomsarray
     */
const assignBonds = function(atoms) {
        // assign bonds - yuck, can't count on connect records

        for (var i = 0, n = atoms.length; i < n; i++) {
            // Don't reindex if atoms are already indexed
            if (!atoms[i].index)
                atoms[i].index = i;
        }

        var grid = {};
        var MAX_BOND_LENGTH = 4.95; // (largest bond length, Cs) 2.25 * 2 * 1.1 (fudge factor)

        for (var index = 0; index < atoms.length; index++) {
            var atom = atoms[index];
            var x = Math.floor(atom.x / MAX_BOND_LENGTH);
            var y = Math.floor(atom.y / MAX_BOND_LENGTH);
            var z = Math.floor(atom.z / MAX_BOND_LENGTH);
            if (!grid[x]) {
                grid[x] = {};
            }
            if (!grid[x][y]) {
                grid[x][y] = {};
            }
            if (!grid[x][y][z]) {
                grid[x][y][z] = [];
            }

            grid[x][y][z].push(atom);
        }

        var findConnections = function(points, otherPoints) {
            for (var i = 0; i < points.length; i++) {
                var atom1 = points[i];
                for (var j = 0; j < otherPoints.length; j++) {
                    var atom2 = otherPoints[j];

                    if (areConnected(atom1, atom2)) {
                        //gracefully handle one-sided bonds
                        var a2i = atom1.bonds.indexOf(atom2.index);
                        var a1i = atom2.bonds.indexOf(atom1.index);
                        if (a2i === -1 && a1i === -1) {
                            atom1.bonds.push(atom2.index);
                            atom1.bondOrder.push(1);
                            atom2.bonds.push(atom1.index);
                            atom2.bondOrder.push(1);
                        } else if (a2i === -1) {
                            atom1.bonds.push(atom2.index);
                            atom1.bondOrder.push(atom2.bondOrder[a1i]);
                        } else if (a1i === -1) {
                            atom2.bonds.push(atom1.index);
                            atom2.bondOrder.push(atom1.bondOrder[a2i]);                          
                        }
                            
                    }
                }
            }
        };


        /*const*/ var OFFSETS = [
            {x: 0, y: 0, z: 1},
            {x: 0, y: 1, z:-1},
            {x: 0, y: 1, z: 0},
            {x: 0, y: 1, z: 1},
            {x: 1, y:-1, z:-1},
            {x: 1, y:-1, z: 0},
            {x: 1, y:-1, z: 1},
            {x: 1, y: 0, z:-1},
            {x: 1, y: 0, z: 0},
            {x: 1, y: 0, z: 1},
            {x: 1, y: 1, z:-1},
            {x: 1, y: 1, z: 0},
            {x: 1, y: 1, z: 1}
        ];
        for (let x in grid) {
            x = parseInt(x);
            for (let y in grid[x]) {
                y = parseInt(y);
                for (let z in grid[x][y]) {
                    z = parseInt(z);
                    let points = grid[x][y][z];

                    for (let i = 0; i < points.length; i++) {
                        let atom1 = points[i];
                        for (let j = i + 1; j < points.length; j++) {
                            let atom2 = points[j];
                            if (areConnected(atom1, atom2)) {
                                if (atom1.bonds.indexOf(atom2.index) === -1) {
                                    atom1.bonds.push(atom2.index);
                                    atom1.bondOrder.push(1);
                                    atom2.bonds.push(atom1.index);
                                    atom2.bondOrder.push(1);
                                }
                            }
                        }
                    }

                    for (let o = 0; o < OFFSETS.length; o++) {
                        let offset = OFFSETS[o];
                        if (!grid[x+offset.x]
                            || !grid[x+offset.x][y+offset.y]
                            || !grid[x+offset.x][y+offset.y][z+offset.z]) continue;

                        let otherPoints = grid[x + offset.x][y + offset.y][z + offset.z];
                        findConnections(points, otherPoints);
                    }
                }
            }
        }
    };

export const chromColors = {

    'chr1' : 0xC0D0FF,        
    'chr2' : 0xB0FFB0,              
    'chr3' : 0xFFC0C8 ,              
    'chr4' : 0xFFFF80,             
    'chr5' : 0xFFC0FF,             
    'chr6' : 0xB0F0F0,            
    'chr7' : 0xFFD070,              
    'chr8' : 0xF08080,        
    'chr9' :  0xF5DEB3,       
    'chr10' :0x00BFFF,             
    'chr11':0xCD5C5C ,              
    'chr12' :0x66CDAA,              
    'chr13' :0x9ACD32,             
    'chr14':0xEE82EE ,            
    'chr15' :  0x00CED1,    
    'chr16': 0x00FF7F,           
    'chr17': 0x3CB371 ,          
    'chr18' :   0x00008B,    
    'chr19' :0xBDB76B,            
    'chr20' :0x006400,            
    'chr21'  :0x800000,     
    'chr22'  :0x808000,       
    'chr23' :0x800080,
    'chrX' :0x008080 ,
     'chrY':0xB8860B,
     'chrM':0xB22222,
    };
/**
 * binary search from https://stackoverflow.com/questions/4431259/formal-way-of-getting-closest-values-in-array-in-javascript-given-a-value-and-a/4431347
 * 
 * @param {*} a: data array
 * @param {*} x: number to be search
 */
export const getClosestValues = (a, x) => {
    if(x < a[0]) {
        return 0;
    }
    if(x > a[a.length-1]) {
        return a.length - 1;
    }
    let lo = -1, hi = a.length;
    while (hi - lo > 1) {
        let mid = Math.round((lo + hi)/2);
        if (a[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    if (a[lo] === x) {
        hi = lo;
    }
    // return [a[lo], a[hi]];
    return lo; // lower index
}