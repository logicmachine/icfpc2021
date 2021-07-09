module.exports = {
  transpileDependencies: [
    'vuetify'
  ],
  publicPath: process.env.NODE_ENV === 'production' ? '/console/' : '/',
  devServer: {
    proxy: 'http://localhost:5000'
  }
}
