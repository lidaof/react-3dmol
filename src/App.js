import React from 'react';
import G3dFile from "g3djs";
import { BigwigSource } from './BigwigSource';
import {chromColors, g3dParser, getClosestValues} from './helpers-3dmol';
import './App.css';



class App extends React.Component {
  constructor(props){
    super(props)
    // this.mol = null;
    // const Parsers = {...window.$3Dmol.Parsers, g3d: func}
    // this.mol = {...window.$3Dmol, Parsers}
    this.mol = window.$3Dmol;
    this.mol.Parsers.g3d = g3dParser;
    this.viewer = null;
    this.viewer2 = null;
    // this.mol.chrom = {};
    // this.mol.chrom.atom = chromColors;
    this.mol.builtinColorSchemes.chrom = {'prop':'chrom', map:chromColors};
    this.myRef = React.createRef();
    this.myRef2 = React.createRef();
    this.state = {
      layout: 'picture',
    };
  }
  
  async componentDidMount(){
    const element = this.myRef.current;
    const element2 = this.myRef2.current;
    const config = { backgroundColor: 'white' };
    this.viewer = this.mol.createViewer( element, {...config, id: 'box1'} );
    this.viewer2 = this.mol.createViewer( element2, {...config, id: 'box2'} );
    this.viewer.linkViewer(this.viewer2);
    this.viewer2.linkViewer(this.viewer);
    const url = "https://target.wustl.edu/dli/tmp/test2.g3d";
    // const url = "https://target.wustl.edu/dli/tmp/k562_1.g3d";
    const file = new G3dFile({ url });
    const data = await file.readData(200000);
    
    console.log(data)
    this.viewer2.addModel( data, "g3d" );
    // this.viewer2.setStyle( {}, { sphere: {colorscheme: 'chrom', opacity: 1, radius: 2}});
    // this.viewer2.setStyle( {}, { stick: {colorscheme: 'chrom', opacity: 1, radius: 0.8}});
    this.viewer2.setStyle({}, { cartoon: {colorscheme: 'chrom', style: 'trace', thickness: 2}});
    // this.viewer2.zoomTo();
    this.viewer2.render();

    this.viewer.addModel( data, "g3d" );
    // this.viewer.setBackgroundColor('white');
    this.viewer.setStyle({}, {line: {colorscheme: 'chrom',opacity: 1}});
    // this.viewer.setViewStyle({style:"outline"});
    // this.viewer.setStyle({chain: 'chr1'}, {cartoon: {colorfunc: colorAsSnake}});
                  // this.viewer.setStyle({},{stick:{colorscheme : 'chain', opacity: 0.4}});
    
    // this.viewer.setStyle( {chain:'chr11'}, { stick: {color: 'grey'}});
    this.viewer.zoomTo();
    this.viewer.render();
    this.viewer.zoom(1.2, 1000);
    
    
    // element.style.border='1px red solid';
    // element2.style.border='1px black solid';
  }


  highlightChrom = () => {
    console.log(this.viewer)
    this.viewer.addStyle( {chain:'chr7'}, { stick: {color: 'grey' }});
  }

  paintWithBigwig = async () => {
    const bw = await this.fetchBwData();
    const bwStarts = bw.map(b=>b.start); //sorted by default
    const scores = bw.map(b=>b.score)
    const maxScore = Math.max(...scores);
    const minScore = Math.min(...scores);
    console.log(maxScore, minScore)
    const grad = new this.mol.Gradient.ROYGB(minScore, maxScore);
    const range = grad.range();
    const colorAsSnake = function(atom) {
      if(atom.properties.start >= 16053398 && atom.properties.start <= 27373766){
        // console.log(atom)
        const a = getClosestValues(bwStarts, atom.properties.start);
        const value = bw[a].score;
        // console.log(a, value)
        return grad.valueToHex(value, range);

      }else {
        return 'grey';
      }
    };
    // console.log(bw, range)
    this.viewer.setStyle( {chain:'chr7'}, { stick: {colorfunc: colorAsSnake }});
  }

  fetchBwData = async () => {
    const bw = new BigwigSource('https://wangftp.wustl.edu/~dli/test/TW463_20-5-bonemarrow_MeDIP.bigWig');
    const bwData = await bw.getData('chr7', 16053398, 27373766, { scale: 1 / 200000 });
    return bwData;
}

onLayoutChange = (e) => {
  this.setState({
    layout: e.target.value
  });
}

  render() {
    return (
      <div className="App">
        <div>
          <p>test 3dmol</p>
          <p><button onClick={this.highlightChrom}>highlight chr7</button></p>
          <p>
            <button onClick={this.paintWithBigwig}>paint with bigwig data ('chr7',26053398,27373766), resolution 200K</button>
          </p>
          <div>
            <strong>Select layout:</strong>
              <ul>
                <li>
                  <label>
                    <input
                      type="radio"
                      value="picture"
                      checked={this.state.layout === "picture"}
                      onChange={this.onLayoutChange}
                    />
                    <span>Picture in picture</span>
                  </label>
                </li>
    
                <li>
                  <label>
                    <input
                      type="radio"
                      value="side"
                      checked={this.state.layout === "side"}
                      onChange={this.onLayoutChange}
                    />
                    <span>Side by side</span>
                  </label>
                </li>
              </ul>
          </div>
        </div>
        <div className={this.state.layout}>
        
        <div className="box1" ref={this.myRef}></div>
        <div className="box2" ref={this.myRef2}></div>
        </div>
      </div>
    );
  }
  
}

export default App;
