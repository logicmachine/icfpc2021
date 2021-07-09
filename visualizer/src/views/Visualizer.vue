<template>
  <v-container fluid fill-height overflow-hidden pa-0>
    <v-row no-gutters>
      <v-col no-butters cols="12" overflow-hidden>
        <svg xmlns="http://www.w3.org/2000/svg"
             :width="canvasWidth" :height="canvasHeight"
             @mousedown.middle="beginScroll"
             @wheel="changeZoom"
             @dragover.prevent
             @drop="dropFile">
          <polygon :points="transformedHolePoints"
                   fill="#bdbdbd" stroke="#000000" stroke-width="2">
          </polygon>
          <line v-for="(edge, index) in figureEdges" :key="'figureEdge' + index"
                :x1="transformX(edge.x0)" :y1="transformY(edge.y0)"
                :x2="transformX(edge.x1)" :y2="transformY(edge.y1)"
                stroke="#616161" stroke-width="2">
          </line>
          <circle v-for="(vertex, index) in figureVertices" :key="'figureVertex' + index"
                  :cx="transformX(vertex.x)" :cy="transformY(vertex.y)" r="4"
                  fill="#616161" stroke="none">
          </circle>
          <line v-for="(edge, index) in solutionEdges" :key="'solutionEdge' + index"
                :x1="transformX(edge.x0)" :y1="transformY(edge.y0)"
                :x2="transformX(edge.x1)" :y2="transformY(edge.y1)"
                stroke="#f44336" stroke-width="2">
          </line>
          <circle v-for="(vertex, index) in solutionVertices" :key="'solutionVertex' + index"
                  :cx="transformX(vertex.x)" :cy="transformY(vertex.y)" r="4"
                  fill="#f44336" stroke="none">
          </circle>
        </svg>
      </v-col>
    </v-row>
  </v-container>
</template>

<script>
export default {
  name: 'Visualizer',

  data: () => {
    return {
      // Environment
      canvasWidth:  0,
      canvasHeight: 0,
      // Viewport
      zoom:    1.0,
      scrollX: 0.0,
      scrollY: 0.0,
      // Scroll
      scrolling:     false,
      beforeScrollX: 0.0,
      beforeScrollY: 0.0,
      scrollOriginX: 0.0,
      scrollOriginY: 0.0,
      // Problem
      problem: {
        "hole": [[45,80],[35,95],[5,95],[35,50],[5,5],[35,5],[95,95],[65,95],[55,80]],
        "epsilon": 150000,
        "figure": {
          "edges": [[2,5],[5,4],[4,1],[1,0],[0,8],[8,3],[3,7],[7,11],[11,13],[13,12],[12,18],[18,19],[19,14],[14,15],[15,17],[17,16],[16,10],[10,6],[6,2],[8,12],[7,9],[9,3],[8,9],[9,12],[13,9],[9,11],[4,8],[12,14],[5,10],[10,15]],
          "vertices":[[20,30],[20,40],[30,95],[40,15],[40,35],[40,65],[40,95],[45,5],[45,25],[50,15],[50,70],[55,5],[55,25],[60,15],[60,35],[60,65],[60,95],[70,95],[80,30],[80,40]]
        }
      },
      // Solution
      solution: {
        "vertices": [[21, 28], [31, 28], [31, 87], [29, 41], [44, 43], [58, 70], [38, 79], [32, 31], [36, 50], [39, 40], [66, 77], [42, 29], [46, 49], [49, 38], [39, 57], [69, 66], [41, 70], [39, 60], [42, 25], [40, 35]]
      }
    }
  },

  created(){
    window.addEventListener('resize', this.resizeHandler)
    window.addEventListener('mouseup', this.mouseupHandler)
    window.addEventListener('mousemove', this.mousemoveHandler)
  },
  destroyed(){
    window.removeEventListener('resize', this.resizeHandler)
    window.removeEventListener('mouseup', this.mouseupHandler)
    window.removeEventListener('mousemove', this.mousemoveHandler)
  },

  mounted(){
    this.resizeHandler()
    this.adjustViewport()
  },

  methods: {
    resizeHandler(){
      this.canvasWidth  = window.innerWidth
      this.canvasHeight = window.innerHeight - 64
    },

    mouseupHandler(e){
      if(e.button == 1){
        this.endScroll(e)
      }
    },

    mousemoveHandler(e){
      this.scroll(e)
    },

    // Drag and drop
    dropFile(e){
      if(e.dataTransfer.files.length != 1){ return }
      const reader = new FileReader()
      reader.onload = e => {
        const obj = JSON.parse(e.target.result)
        const hasProperty = name => { return Object.prototype.hasOwnProperty.call(obj, name) }
        if(hasProperty('hole') && hasProperty('figure')){
          this.solution = { vertices: [] }
          this.problem  = obj
        }else if(hasProperty('vertices')){
          if(this.problem.figure.vertices.length == obj.vertices.length){
            this.solution = obj
          }
        }
      }
      reader.readAsText(e.dataTransfer.files[0])
      e.preventDefault()
    },

    // Scroll
    beginScroll(e){
      this.scrolling = true
      this.beforeScrollX = this.scrollX
      this.beforeScrollY = this.scrollY
      this.scrollOriginX = e.screenX
      this.scrollOriginY = e.screenY
      e.preventDefault()
    },

    scroll(e){
      if(this.scrolling){
        const deltaX = e.screenX - this.scrollOriginX
        const deltaY = e.screenY - this.scrollOriginY
        this.scrollX = this.beforeScrollX - deltaX
        this.scrollY = this.beforeScrollY - deltaY
      }
    },

    endScroll(){
      this.scrolling = false
    },

    // Zoom
    changeZoom(e){
      const scale = e.deltaY > 0 ? 0.5 : 2.0
      const offsetX = e.offsetX, offsetY = e.offsetY
      const imageX = offsetX + this.scrollX
      const imageY = offsetY + this.scrollY
      this.scrollX = imageX * scale - offsetX
      this.scrollY = imageY * scale - offsetY
      this.zoom *= scale
      e.preventDefault()
    },

    // Viewport transform
    transformX(x){ return x * this.zoom - this.scrollX },
    transformY(y){ return y * this.zoom - this.scrollY },
    transformSize(s){ return s * this.zoom },

    // Adjust viewport: adjust zoom and centering
    adjustViewport(){
      var minx = Infinity, maxx = -Infinity
      var miny = Infinity, maxy = -Infinity
      const accumulate = a => {
        minx = Math.min(minx, a[0])
        maxx = Math.max(maxx, a[0])
        miny = Math.min(miny, a[1])
        maxy = Math.max(maxy, a[1])
      }
      this.problem.hole.forEach(accumulate)
      this.problem.figure.vertices.forEach(accumulate)
      if(minx > maxx || miny > maxy){ return }
      const w = maxx - minx
      const h = maxy - miny
      const z = Math.min(this.canvasWidth / w, this.canvasHeight / h) * 0.9
      this.scrollX = (maxx + minx) * 0.5 * z - this.canvasWidth  * 0.5
      this.scrollY = (maxy + miny) * 0.5 * z - this.canvasHeight * 0.5
      this.zoom    = z
    },
  },

  computed: {
    hole(){
      return this.problem.hole.map(v => {
        return { x: v[0], y: v[1] }
      })
    },
    transformedHolePoints(){
      const points = []
      this.problem.hole.forEach(v => {
        const x = this.transformX(v[0])
        const y = this.transformY(v[1])
        points.push(`${x},${y}`)
      })
      return points.join(' ')
    },
    figureEdges(){
      const vertices = this.problem.figure.vertices
      return this.problem.figure.edges.map(e => {
        const u = vertices[e[0]], v = vertices[e[1]]
        return { x0: u[0], y0: u[1], x1: v[0], y1: v[1] }
      })
    },
    figureVertices(){
      return this.problem.figure.vertices.map(v => {
        return { x: v[0], y: v[1] }
      })
    },
    solutionEdges(){
      const vertices = this.solution.vertices
      if(!vertices || vertices.length == 0){ return [] }
      return this.problem.figure.edges.map(e => {
        const u = vertices[e[0]], v = vertices[e[1]]
        return { x0: u[0], y0: u[1], x1: v[0], y1: v[1] }
      })
    },
    solutionVertices(){
      return this.solution.vertices.map(v => {
        return { x: v[0], y: v[1] }
      })
    },
  }
}
</script>

<style>
html {
  overflow: hidden;
}
</style>
