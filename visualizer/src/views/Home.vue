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
          this.problems = response.data
        })
    }
  }
}
</script>
