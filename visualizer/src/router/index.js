import Vue from 'vue'
import VueRouter from 'vue-router'
import Home from '../views/Home.vue'
import Visualizer from '../views/Visualizer.vue'

Vue.use(VueRouter)

const routes = [
  {
    path: '/',
    name: 'Home',
    component: Home
  },
  {
    path: '/visualize',
    name: 'Visualizer',
    component: Visualizer
  },
  {
    path: '/visualize/:id',
    name: 'Visualizer',
    component: Visualizer
  }
]

const router = new VueRouter({
  mode: 'hash',
  base: process.env.BASE_URL,
  routes
})

export default router
