<template>
  <v-container>
    <v-row>
      <v-col>
        <h2 class="headline mb-3">Problems</h2>
        <v-data-table :headers="headers" :items="problems" class="elevation-1" :options="options">
          <template #item.id="{ item }">
            <router-link :to="{ name: 'SubmissionList', params: { problem_id: item.id } }">
              {{ item.id }}
            </router-link>
          </template>
          <template #item.name="{ item }">
            {{ item.name }}
            <v-chip :color="bonusColor(item.bonus)" small class="ma-2">{{ bonusLabel(item.bonus) }}</v-chip>
          </template>
          <template #item.best_dislikes="{ item }">
            {{ (item.best_dislikes == -1) ? "❌" : item.best_dislikes}}
          </template>
        </v-data-table>
      </v-col>
    </v-row>
  </v-container>
</template>

<script>
import axios from 'axios'

const api_options = {
  headers: {
    'Authorization': `Bearer ${process.env.VUE_APP_API_TOKEN}`
  }
}

export default {
  name: 'Home',

  data: () => {
    return {
      problem_id: null,
      headers: [
        { text: 'ID', value: 'id' },
        { text: 'Name', value: 'name' },
        { text: 'Best', value: 'best_dislikes' },
        { text: 'Min', value: 'min_dislikes' },
        { text: 'Our Score', value: 'score' },
        { text: 'Max Score', value: 'max_score' },
        { text: '|hole|', value: 'hole_size' },
        { text: '|edge|', value: 'num_edges' },
        { text: '|vertices|', value: 'num_vertices' },
      ],
      problems: [],
      options: {
        itemsPerPage: -1
      }
    }
  },

  created(){
    this.fetchData()
  },

  watch: {
    $route: 'fetchData'
  },

  methods: {
    fetchData(){
      axios
        .get(`/api/problems`, api_options)
        .then(response => {
          response.data.forEach(item => {
            const limit = 1000 * Math.log2(item.hole_size * item.num_edges * item.num_vertices / 6)
            if(item.best_dislikes === null || item.best_dislikes >= 9e18){ item.best_dislikes = -1 }
            item.max_score = Math.ceil(limit)
            if(item.best_dislikes < 0){
              item.score = 0
            }else{
              item.score = Math.ceil(limit * Math.sqrt((item.min_dislikes + 1) / (item.best_dislikes + 1)))
            }
          })
          this.problems = response.data
        })
    },

    bonusLabel(type){
      if(type == 0){
        return 'None'
      }else if(type == 1){
        return 'Globalist'
      }else if(type == 2){
        return 'Break a Leg'
      }else if(type == 3){
        return 'Wallhack'
      }else if(type == 4){
        return 'Superflex'
      }
    },

    bonusColor(type){
      if(type == 0){
        return 'gray'
      }else if(type == 1){
        return 'yellow'
      }else if(type == 2){
        return 'blue'
      }else if(type == 3){
        return 'orange'
      }else if(type == 4){
        return 'cyan'
      }
    }
  }
}
</script>
